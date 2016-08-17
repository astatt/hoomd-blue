// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

#include "hoomd/HOOMDMath.h"
#include "MuellerPlatheFlow.h"
#include "hoomd/HOOMDMPI.h"

namespace py=pybind11;
using namespace std;

const unsigned int INVALID_TAG=UINT_MAX;
const Scalar INVALID_VEL = FLT_MAX; //should be ok, even for double.


//! \file MuellerPlatheFlow.cc Implementation of CPU version of MuellerPlatheFlow.

MuellerPlatheFlow::MuellerPlatheFlow(std::shared_ptr<SystemDefinition> sysdef,
                                     std::shared_ptr<ParticleGroup> group,
                                     std::shared_ptr<Variant> flow_target,
                                     const unsigned int slab_direction,
                                     const unsigned int flow_direction,
                                     const unsigned int N_slabs,
                                     const unsigned int min_slab,
                                     const unsigned int max_slab)
    :Updater(sysdef), m_group(group), m_flow_target(flow_target)
    ,m_flow_epsilon(1e-2)
    ,m_N_slabs(N_slabs),m_min_slab(min_slab),m_max_slab(max_slab)
    ,m_exchanged_momentum(0),m_has_min_slab(true),m_has_max_slab(true)
    {
    assert(m_flow_target);
    switch(slab_direction)
        {
        case 0: m_slab_direction = X; break;
        case 1: m_slab_direction = Y; break;
        case 2: m_slab_direction = Z; break;
        default:
        throw runtime_error("ERROR: invalid slab direction.");
        }
    switch(flow_direction)
        {
        case 0: m_flow_direction = X; break;
        case 1: m_flow_direction = Y; break;
        case 2: m_flow_direction = Z; break;
        default:
        throw runtime_error("ERROR: invalid slow direction.");
        }

    m_last_max_vel.s = -INVALID_VEL;
    m_last_max_vel.i = INVALID_TAG;
    m_last_min_vel.s = INVALID_VEL;
    m_last_min_vel.i = INVALID_TAG;

    m_exec_conf->msg->notice(5) << "Constructing MuellerPlatheFlow " << endl;
    this->update_domain_decomposition();
    //m_exec_conf->msg->notice(0)<<m_exec_conf->getRank()<<": "<< m_max_swap.gbl_rank<<" "<<m_min_swap.gbl_rank<<endl;

    //Check min max slab.
    this->set_min_slab(m_min_slab);
    this->set_max_slab(m_max_slab);
    }

MuellerPlatheFlow::~MuellerPlatheFlow(void)
    {
    m_exec_conf->msg->notice(5) << "Destroying MuellerPlatheFlow " << endl;
    }

void MuellerPlatheFlow::update(unsigned int timestep)
    {
    const BoxDim&box= m_pdata->getBox();
    double area;
    switch(m_slab_direction)
      {
      case X: area = box.getL().y * box.getL().z;
      case Y: area = box.getL().x * box.getL().z;
      case Z: area = box.getL().y * box.getL().x;
      }
    unsigned int counter = 0;
    const unsigned int max_iteration = 100;
    while( fabs( m_flow_target->getValue(timestep) -
                 this->summed_exchanged_momentum()/area ) > this->get_flow_epsilon()
           && counter < max_iteration)
        {
        counter++;
        if( m_flow_target->getValue(timestep)*
            this->summed_exchanged_momentum() < 0)
            {
            this->swap_min_max_slab();
            }

        m_last_max_vel.s = -INVALID_VEL;
        m_last_max_vel.i = INVALID_TAG;
        m_last_min_vel.s = INVALID_VEL;
        m_last_max_vel.i = INVALID_TAG;
        search_min_max_velocity();
#ifdef ENABLE_MPI
        mpi_exchange_velocity();

#endif//ENABLE_MPI
        if( m_last_max_vel.s == -INVALID_VEL
            || static_cast<unsigned int>(m_last_max_vel.i) == INVALID_TAG
            || m_last_min_vel.s == -INVALID_VEL
            || static_cast<unsigned int>(m_last_min_vel.i) == INVALID_TAG)
            {
            m_exec_conf->msg->warning() << "WARNING: at time "<<timestep
                                            <<"  MuellerPlatheFlow could not find a min/max pair."
                                        <<" This could cause infinite loops."<<endl;
            }
        else
            {
            update_min_max_velocity();
            const int sign = this->get_max_slab() > this->get_min_slab() ? 1 : -1;
            m_exchanged_momentum += sign * (m_last_max_vel.s - m_last_min_vel.s);
            }
        }
    if(counter >= max_iteration)
      {
        stringstream s;
        s<<" After "<<counter<<" MuellerPlatheFlow could not achieve the target: "
         <<m_flow_target->getValue(timestep)<<" only "
         <<this->summed_exchanged_momentum()/area<<" could be achieved."<<endl;
        m_exec_conf->msg->warning()<<s.str();
      }
    // stringstream s;
    // s<<this->summed_exchanged_momentum()/area<<"\t"<<m_flow_target->getValue(timestep)<<endl;
    // m_exec_conf->msg->collectiveNoticeStr(0,s.str());
    }

void MuellerPlatheFlow::swap_min_max_slab(void)
    {
    std::swap( m_max_slab, m_min_slab);
//     const unsigned int old_max_slab = this->get_max_slab();
//     m_max_slab = this->get_min_slab();
//     m_min_slab = old_max_slab;

    std::swap( m_last_max_vel, m_last_min_vel);
//     const Scalar_Int old_max_vel = m_last_max_vel;
//     m_last_max_vel = m_last_min_vel;
//     m_last_min_vel = old_max_vel;

    std::swap( m_has_max_slab, m_has_min_slab);
//     const bool old_max_bool = m_has_max_slab;
//     m_has_max_slab = m_has_min_slab;
//     m_has_min_slab = old_max_bool;

#ifdef ENABLE_MPI
    std::swap( m_max_swap, m_min_swap);
//     const struct MPI_SWAP old_max_swap = m_max_swap;
//     m_max_swap = m_min_swap;
//     m_min_swap = old_max_swap;
#endif//ENABLE_MPI
    }

void MuellerPlatheFlow::set_min_slab(const unsigned int min_slab)
    {
    if(min_slab >= m_N_slabs)
        {
        m_exec_conf->msg->error()<< "MuellerPlatheFlow is initialized with invalid min_slab: "
                        <<min_slab<<"/"<<m_N_slabs<<endl;
        throw runtime_error("ERROR: Invalid min_slab.\n");
        }
    if( min_slab != m_min_slab)
        this->update_domain_decomposition();
    }

void MuellerPlatheFlow::set_max_slab(const unsigned int max_slab)
    {
    if(max_slab >= m_N_slabs)
        {
        m_exec_conf->msg->error()<< "MuellerPlatheFlow is initialized with invalid max_slab: "
                        <<max_slab<<"/"<<m_N_slabs<<endl;
        throw runtime_error("ERROR: Invalid max_slab.\n");
        }
    if(max_slab != m_max_slab)
        this->update_domain_decomposition();
    }

void MuellerPlatheFlow::update_domain_decomposition(void)
    {
#ifdef ENABLE_MPI
    std::shared_ptr<DomainDecomposition> dec = m_pdata->getDomainDecomposition();
    if( dec )
        {
        const Scalar min_frac = m_min_slab/static_cast<Scalar>(m_N_slabs);
        const Scalar max_frac = m_max_slab/static_cast<Scalar>(m_N_slabs);

        const uint3 grid = dec->getGridPartition();
        const uint3 pos = dec->getGridPos();
        const unsigned int my_grid = m_slab_direction == X ? grid.x :( m_slab_direction == Y ? grid.y : grid.z);
        const unsigned int my_pos = m_slab_direction == X ? pos.x :( m_slab_direction == Y ? pos.y : pos.z);
        if( m_N_slabs % my_grid != 0)
            {
            m_exec_conf->msg->warning()<<"MuellerPlatheFlow::N_slabs does is divideable "
                " by the domain decomposition. Adjusting N_slabs.\n"<<endl;
            m_N_slabs += m_N_slabs % my_grid;
            m_min_slab = static_cast<unsigned int>(min_frac*m_N_slabs);
            m_max_slab = static_cast<unsigned int>(max_frac*m_N_slabs);
            }

        m_has_min_slab = false;
        if( my_pos == this->get_min_slab() / (m_N_slabs/my_grid) )
            m_has_min_slab = true;

        m_has_max_slab = false;
        if( my_pos == this->get_max_slab() / (m_N_slabs/my_grid) )
            m_has_max_slab = true;


        //Create the communicator.
        const int min_color =  this->has_min_slab() ? 0 : MPI_UNDEFINED;
        init_mpi_swap(&m_min_swap,min_color);

        const int max_color =  this->has_max_slab() ? 0 : MPI_UNDEFINED;
        init_mpi_swap(&m_max_swap,max_color);
        }
#endif//ENABLE_MPI
    }

void MuellerPlatheFlow::search_min_max_velocity(void)
    {
    const unsigned int group_size = m_group->getNumMembers();
    if(group_size == 0)
        return;
    if( !this->has_max_slab() and !this->has_min_slab())
        return;
    if(m_prof) m_prof->push("MuellerPlatheFlow::search");
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(),access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(),access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(),access_location::host, access_mode::read);
    const BoxDim& gl_box = m_pdata->getGlobalBox();
    for(unsigned int group_idx=0; group_idx < group_size; group_idx++)
        {
        const unsigned int j = m_group->getMemberIndex(group_idx);
        if( j < m_pdata->getN() )
            {
            unsigned int index;
            switch(m_slab_direction)
                {
                case X: index = (( (h_pos.data[j].x)/gl_box.getL().x + .5) * this->get_N_slabs());
                    break;
                case Y: index = ( (h_pos.data[j].y)/gl_box.getL().y + .5) * this->get_N_slabs();
                    break;
                case Z: index = ( (h_pos.data[j].z)/gl_box.getL().z + .5) * this->get_N_slabs();
                    break;
                }
            index %= this->get_N_slabs(); //border cases. wrap periodic box
            assert( index >= 0);
            assert( index < this->get_N_slabs());
            if( index == this->get_max_slab() || index == this->get_min_slab())
                {
                Scalar vel;
                switch(m_flow_direction)
                    {
                    case X: vel = h_vel.data[j].x;break;
                    case Y: vel = h_vel.data[j].y;break;
                    case Z: vel = h_vel.data[j].z;break;
                    }
                if( index == this->get_max_slab() && m_last_max_vel.s < vel &&  this->has_max_slab())
                    {
                    m_last_max_vel.s = vel;
                    m_last_max_vel.i = h_tag.data[j];
                    }
                if( index == this->get_min_slab() && m_last_min_vel.s > vel && this->has_min_slab())
                    {
                    m_last_min_vel.s = vel;
                    m_last_min_vel.i = h_tag.data[j];
                    }
                }
            }
        }
    if(m_prof) m_prof->pop();

    }

void MuellerPlatheFlow::update_min_max_velocity(void)
    {
    if(m_prof) m_prof->push("MuellerPlatheFlow::update");
    ArrayHandle<unsigned int> h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);
    const unsigned int min_idx = h_rtag.data[m_last_min_vel.i];
    const unsigned int max_idx = h_rtag.data[m_last_max_vel.i];
    const unsigned int Ntotal = m_pdata->getN()+m_pdata->getNGhosts();
    //Is my particle local on the processor?
    if( min_idx < Ntotal || max_idx < Ntotal)
        {
        ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(),access_location::host,access_mode::readwrite);
        //Swap the particles the new velocities.
        if( min_idx < Ntotal)
            {
            switch(m_flow_direction)
                {
                case X: h_vel.data[min_idx].x = m_last_max_vel.s; break;
                case Y: h_vel.data[min_idx].y = m_last_max_vel.s; break;
                case Z: h_vel.data[min_idx].z = m_last_max_vel.s; break;
                }
            }
        if( max_idx < Ntotal)
            switch(m_flow_direction)
                {
                case X: h_vel.data[max_idx].x = m_last_min_vel.s;
                case Y: h_vel.data[max_idx].y = m_last_min_vel.s;
                case Z: h_vel.data[max_idx].z = m_last_min_vel.s;
                }
        }
    if(m_prof) m_prof->pop();
    }

#ifdef ENABLE_MPI

//Not performance optimized: could be slow. It is meant for init.
void MuellerPlatheFlow::init_mpi_swap(struct MPI_SWAP*ms,const int color)
    {
    if( ms->initialized && ms->rank != MPI_UNDEFINED)
        MPI_Comm_free( &(ms->comm) );

    MPI_Comm_split( m_exec_conf->getMPICommunicator(),
                    color, m_exec_conf->getRank(),
                    &(ms->comm));
    if( color != MPI_UNDEFINED)
        {
        MPI_Comm_rank(ms->comm,&(ms->rank) );
        MPI_Comm_size(ms->comm,&(ms->size) );
        }
    else
        {
        ms->rank = MPI_UNDEFINED;
        ms->size = MPI_UNDEFINED;
        }
    int send =0;
    if( ms->rank == 0)
        {
        send = 1;
        }
    vector<int> recv;
    all_gather_v(send,recv,m_exec_conf->getMPICommunicator());

    //ARG_MAX;
    ms->gbl_rank = -1;
    for(unsigned int i=0; i < recv.size(); i++)
        {
        if( recv[i] == 1)
            {
            assert( ms->gbl_rank == -1);
            ms->gbl_rank = i;
            }
        }
    assert( ms->gbl_rank > -1);
    ms->initialized = true;
    }

void MuellerPlatheFlow::bcast_vel_to_all(struct MPI_SWAP*ms,Scalar_Int*vel,const MPI_Op op)
    {
    if( ms->rank != MPI_UNDEFINED)
        {
        Scalar_Int tmp;
        tmp.s = vel->s;
        tmp.i = ms->rank;
        MPI_Allreduce(MPI_IN_PLACE,&tmp,1,MPI_HOOMD_SCALAR_INT,op, ms->comm);
        //Send to subset local rank
        if( tmp.i != 0)
            {
            if( ms->rank == 0)
                {
                MPI_Recv( vel,1, MPI_HOOMD_SCALAR_INT, tmp.i, 42, ms->comm, MPI_STATUS_IGNORE);
                }
            if( ms->rank == tmp.i)
                {
                MPI_Send( vel,1, MPI_HOOMD_SCALAR_INT, 0, 42, ms->comm);
                }
            }
        //ms->rank == 0 has the definite answer
        }
    //Broadcast the result to every rank.
    //This way, each rank can check, whether
    //it needs to update and can determine the exchanged momentum.
    MPI_Bcast( vel, 1, MPI_HOOMD_SCALAR_INT, ms->gbl_rank, m_exec_conf->getMPICommunicator());
    }

void MuellerPlatheFlow::mpi_exchange_velocity(void)
    {
#ifdef ENABLE_MPI
    if( m_pdata->getDomainDecomposition() )
        {
        bcast_vel_to_all(&m_min_swap,&m_last_min_vel,MPI_MINLOC);
        bcast_vel_to_all(&m_max_swap,&m_last_max_vel,MPI_MAXLOC);
        }
#endif//ENABLE_MPI
    }

#endif//ENABLE_MPI

void export_MuellerPlatheFlow(py::module& m)
    {
    py::class_< MuellerPlatheFlow, std::shared_ptr<MuellerPlatheFlow> >
        (m,"MuellerPlatheFlow",py::base<Updater>())
        .def(py::init< std::shared_ptr<SystemDefinition>,std::shared_ptr<ParticleGroup>,
             std::shared_ptr<Variant>, const unsigned int, const unsigned int,
             const unsigned int, const unsigned int, const unsigned int >() )
        .def("getNSlabs",&MuellerPlatheFlow::get_N_slabs)
        .def("getMinSlab",&MuellerPlatheFlow::get_min_slab)
        .def("getMaxSlab",&MuellerPlatheFlow::get_max_slab)
        .def("hasMinSlab",&MuellerPlatheFlow::has_min_slab)
        .def("hasMaxSlab",&MuellerPlatheFlow::has_max_slab)
        .def("getFlowEpsilon",&MuellerPlatheFlow::get_flow_epsilon)
        .def("setFlowEpsilon",&MuellerPlatheFlow::set_flow_epsilon)
        .def("getSummedExchangedMomentum",&MuellerPlatheFlow::summed_exchanged_momentum)
        // Functions not needed for pyhton interface users.
        // .def("setMinSlab",&MuellerPlatheFlow::set_min_slab)
        // .def("setMaxSlab",&MuellerPlatheFlow::set_max_slab)
        // .def("swapMinMaxSlab",&MuellerPlatheFlow::swap_min_max_slab)
        // .def("updateDomainDecomposition",&MuellerPlatheFlow::update_domain_decomposition)
        ;
    }
