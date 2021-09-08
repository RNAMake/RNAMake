////
////  resource_manager.cc
////  RNAMake
////
////  Created by Joseph Yesselman on 8/8/15.
////  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
////
//
//#include "structure/residue_type_set_manager.h"
//#include "resources/resource_manager.h"
//
//namespace resources {
//
//Manager::Manager() {
//    mf_ = motif::MotifFactory();
//    added_motifs_ = AddedMotifLibrary();
//    mlibs_ = std::map<String, MotifSqliteLibraryOP>();
//    ms_libs_ = std::map<String, MotifStateSqliteLibraryOP>();
//    mse_libs_ = std::map<String, MotifStateEnsembleSqliteLibraryOP>();
//    extra_me_ = std::map<String, motif::MotifEnsembleOP>();
//
//    for (auto const & kv : MotifSqliteLibrary::get_libnames()) {
//        //if(kv.first == "bp_steps") { continue; }
//        mlibs_[kv.first] = std::make_shared<MotifSqliteLibrary>(kv.first);
//    }
//
//    for (auto const & kv : MotifStateSqliteLibrary::get_libnames()) {
//        //if(kv.first == "bp_steps") { continue; }
//        ms_libs_[kv.first] = std::make_shared<MotifStateSqliteLibrary>(kv.first);
//    }
//
//    for (auto const & kv : MotifStateEnsembleSqliteLibrary::get_libnames()) {
//        mse_libs_[kv.first] = std::make_shared<MotifStateEnsembleSqliteLibrary>(kv.first);
//    }
//
//}
//
//
//// getting functions  //////////////////////////////////////////////////////////////////////////////
//
//motif::MotifOP
//Manager::bp_step(
//        String const & end_id) {
//    return mlibs_["new_bp_steps"]->get("", end_id);
//}
//
//motif::MotifStateOP
//Manager::bp_step_state(
//        String const & end_id) {
//    return bp_step(end_id)->get_state();
//}
//
//motif::MotifOP
//Manager::motif(
//        String const & name,
//        String const & end_id,
//        String const & end_name) {
//
//    for (auto const & kv : mlibs_) {
//        if (kv.second->contains(name, end_id, end_name)) {
//            return kv.second->get(name, end_id, end_name);
//        }
//    }
//
//    if (added_motifs_.contains(name, end_id, end_name)) {
//        return added_motifs_.get(name, end_id, end_name);
//    }
//
//    throw ResourceManagerException(
//            "cannot find motif in resource manager with search: "
//                    "name=" + name + " end_id=" + end_id + " end_name= " + end_name);
//
//}
//
//motif::MotifStateOP
//Manager::motif_state(
//        String const & name,
//        String const & end_id,
//        String const & end_name) {
//
//    for (auto const & kv : ms_libs_) {
//        if (kv.second->contains(name, end_id, end_name)) {
//            return kv.second->get(name, end_id, end_name);
//        }
//    }
//
//    if (added_motifs_.contains(name, end_id, end_name)) {
//        return added_motifs_.get(name, end_id, end_name)->get_state();
//    }
//
//    throw ResourceManagerException(
//            "cannot find motif state in resource manager with search: "
//                    "name=" + name + " end_id=" + end_id + " end_name= " + end_name);
//}
//
//motif::MotifStateEnsembleOP
//Manager::motif_state_ensemble(
//        String const & name) {
//
//    for (auto const & kv : mse_libs_) {
//        if (kv.second->contains(name)) {
//            return kv.second->get(name);
//        }
//    }
//
//    throw ResourceManagerException(
//            "cannot find motif state ensemble in resource manager with search: "
//                    "name=" + name);
//}
//
//
//// add functions  //////////////////////////////////////////////////////////////////////////////////
//
//structure::RNAStructureOP
//Manager::get_structure(
//        String const & path,
//        String name,
//        int force_num_chains) {
//
//    auto m = mf_.motif_from_file(path, false, false, force_num_chains);
//    if (name != "") { m->name(name); }
//
//    for (auto & end : m->ends()) {
//        int flip_res = 0;
//        for (auto const & c : m->chains()) {
//            if (c->first() == end->res2()) {
//                flip_res = 1;
//                break;
//            }
//        }
//
//        if (!flip_res) { continue; }
//
//        auto temp = end->res1();
//        end->res1(end->res2());
//        end->res2(temp);
//    }
//
//    return m;
//}
//
//
//void
//Manager::add_motif(
//        String const & path,
//        String name,
//        util::MotifType mtype) {
//
//    auto m = mf_.motif_from_file(path, false, true);
//    m->mtype(mtype);
//
//    if (name != "") { m->name(name); }
//
//    motif::MotifOPs motifs;
//    std::map<util::Uuid, String, util::UuidCompare> end_ids;
//    for (int i = 0; i < m->ends().size(); i++) {
//        auto m_added = mf_.get_oriented_motif(m, i);
//
//        motifs.push_back(m_added);
//        end_ids[m_added->ends()[0]->uuid()] = m_added->end_ids()[0];
//    }
//
//    if (motifs.size() == 0) {
//        throw ResourceManagerException(
//                "attempted to add motif from path " + path + " unforunately it has no viable "
//                        "basepair ends to be build from ");
//    }
//
//    for (auto & m : motifs) {
//        Strings final_end_ids;
//        for (auto const & end : m->ends()) {
//            final_end_ids.push_back(end_ids[end->uuid()]);
//        }
//        m->end_ids(final_end_ids);
//        added_motifs_.add_motif(m);
//    }
//
//}
//
//void
//Manager::add_motif(
//        motif::MotifOP const & m,
//        String name) {
//
//    if (name != "") { m->name(name); }
//
//    auto motifs = motif::MotifOPs();
//    auto end_ids = std::map<util::Uuid, String, util::UuidCompare>();
//    for (int i = 0; i < m->ends().size(); i++) {
//        auto m_added = mf_.can_align_motif_to_end(m, i);
//        if (m_added == nullptr) {
//            //TODO switch this to logging system
//            std::cout << "RESOURCE MANAGER WARNING: cannot create standardized motif for ";
//            std::cout << m->name() << " with end" << m->ends()[i]->name() << std::endl;
//            continue;
//        }
//        m_added = mf_.align_motif_to_common_frame(m_added, i);
//        if (m_added == nullptr) {
//            std::cout << "RESOURCE MANAGER WARNING: cannot create standardized motif for ";
//            std::cout << m->name() << " with end" << m->ends()[i]->name() << std::endl;
//            continue;
//        }
//
//        motifs.push_back(m_added);
//        end_ids[m_added->ends()[0]->uuid()] = m_added->end_ids()[0];
//    }
//
//    if (motifs.size() == 0) {
//        throw ResourceManagerException(
//                "attempted to add motif " + m->name() + " unforunately it has no viable "
//                        "basepair ends to be build from ");
//    }
//
//    for (auto & m : motifs) {
//        auto final_end_ids = Strings();
//        for (auto const & end : m->ends()) {
//            final_end_ids.push_back(end_ids[end->uuid()]);
//        }
//        m->end_ids(final_end_ids);
//        added_motifs_.add_motif(m);
//    }
//}
//
//
//void
//Manager::register_motif(
//        motif::MotifOP const & m) {
//
//    if(m->name() == "") {
//        throw ResourceManagerException(
//            "attempted to register motif with no name this will make it "
//            "extremely unlikely you will be able to retrieve it properly!");
//    }
//
//    added_motifs_.add_motif(m);
//
//
//}
//
//void
//Manager::register_extra_motif_ensembles(
//        String const & f_name) {
//
//    auto lines = base::get_lines_from_file(f_name);
//
//    for (auto const & l : lines) {
//        if (l.length() < 10) { continue; }
//        auto spl = base::split_str_by_delimiter(l, "!!");
//        extra_me_[spl[0]] = std::make_shared<motif::MotifEnsemble>(spl[1],
//                                                                   structure::ResidueTypeSetManager::getInstance().residue_type_set());
//
//        for (auto const & mem : extra_me_[spl[0]]->members()) {
//            try {
//                added_motifs_.add_motif(mem->motif);
//            } catch (AddedMotifLibraryException) {}
//        }
//
//        std::cout << "RESOURCE MANAGER: motif ensemble for " << spl[0] << " registered" << std::endl;
//    }
//
//}
//
//int
//Manager::has_supplied_motif_ensemble(
//        String const & m_name,
//        String const & end_name) {
//
//    auto key = m_name + "-" + end_name;
//    if (extra_me_.find(key) == extra_me_.end()) { return 0; }
//    else { return 1; }
//
//}
//
//motif::MotifEnsembleOP const &
//Manager::get_supplied_motif_ensemble(
//        String const & m_name,
//        String const & end_name) {
//
//    auto key = m_name + "-" + end_name;
//    if (extra_me_.find(key) == extra_me_.end()) {
//        throw ResourceManagerException(
//                "motif ensemble with name: " + m_name + " and end_name: " + end_name +
//                " has not been supplied please use register_extra_motif_ensembles to "
//                        "do so");
//    }
//
//    return extra_me_[key];
//
//}
//
//}
//
//
//
//
//
//
//
//
//
//
//
//
