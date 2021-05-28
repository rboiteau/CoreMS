__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

from copy import deepcopy
import itertools
import multiprocessing
import json
import cProfile
import io
import pstats
import contextlib

from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import load_only
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, func
from tqdm import tqdm

from corems.encapsulation.factory.processingSetting import MolecularLookupDictSettings
from corems.encapsulation.constant import Atoms
from corems.molecular_id.factory.molecularSQL import CarbonHydrogen, HeteroAtoms, MolecularFormulaLink
from corems.encapsulation.factory.parameters import MSParameters
from corems import chunks, timeit
from corems.molecular_id.factory.molecularSQL import MolForm_SQL
import os

@contextlib.contextmanager
def profiled():
    pr = cProfile.Profile()
    pr.enable()
    yield
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    ps.print_stats()
    # uncomment this to see who's calling what
    # ps.print_callers()
    print(s.getvalue())

def insert_database_worker(args):
        
        results, url = args
        
        if not url:
            
            url = 'sqlite:///db/molformulas.sqlite'

        if url[0:6] == 'sqlite':
            engine = create_engine(url, echo = False)
        else:
            engine = create_engine(url, echo = False, isolation_level="AUTOCOMMIT")
        
        session_factory = sessionmaker(bind=engine)
        session = session_factory()
        insert_query = MolecularFormulaLink.__table__.insert().values(results)
        session.execute(insert_query)
        session.commit()
        session.close()
        engine.dispose()

class MolecularCombinations:
     
    def __init__(self, sql_db = None):

        if not sql_db:
            
            self.sql_db = MolForm_SQL()
        else:
            
            self.sql_db = sql_db

    def cProfile_worker(self, args):
        
        cProfile.runctx('self.get_mol_formulas(*args)', globals(), locals(), 'mf_database_cprofile.prof')

    def check_database_get_class_list(self, molecular_search_settings):
        
        all_class_to_create = []
        
        classes_dict = self.get_classes_in_order(molecular_search_settings)
        
        class_str_set = set(classes_dict.keys())
        
        existing_classes_objs = self.sql_db.session.query(HeteroAtoms).distinct().all()
        
        existing_classes_str = set([classe.name for classe in existing_classes_objs])

        self.len_existing_classes = len(existing_classes_str)

        class_to_create = class_str_set - existing_classes_str
        
        class_count= len(existing_classes_objs)

        data_classes = list()    
        for index, class_str in enumerate(class_to_create):
            
            class_dict = classes_dict.get(class_str)
            halogen_count = self.get_total_halogen_atoms(class_dict)
            data_classes.append({"name":class_str, "id":class_count+ index + 1, "halogensCount": halogen_count})

        #data_classes = [{"name":class_str, "id":class_count+ index + 1} for index, class_str in enumerate(class_to_create)]

        if data_classes:

            list_insert_chunks = chunks(data_classes, self.sql_db.chunks_count)
            for insert_chunk in  list_insert_chunks:   
                insert_query = HeteroAtoms.__table__.insert().values(insert_chunk)
                self.sql_db.session.execute(insert_query)

        for index, class_str in enumerate(class_to_create):

            class_tuple = (class_str, classes_dict.get(class_str), class_count+ index + 1)

            all_class_to_create.append(class_tuple)

        return [(c_s, c_d) for c_s, c_d in classes_dict.items()], all_class_to_create, existing_classes_objs       

    def get_carbonsHydrogens(self, settings, odd_even):

        operator = '==' if odd_even == 'even' else '!=' 
        usedAtoms = settings.usedAtoms
        user_min_c, user_max_c = usedAtoms.get('C')
        user_min_h, user_max_h = usedAtoms.get('H')

        return eval("self.sql_db.session.query(CarbonHydrogen).filter(" 
                                       "CarbonHydrogen.C >= user_min_c,"
                                        "CarbonHydrogen.H >= user_min_h,"
                                        "CarbonHydrogen.C <= user_max_c,"
                                        "CarbonHydrogen.H <= user_max_h,"
                                        "CarbonHydrogen.H % 2" + operator+ "0).all()")

    def add_carbonsHydrogens(self, settings, existing_classes_objs):

        usedAtoms = settings.usedAtoms

        user_min_c, user_max_c = usedAtoms.get('C')
        user_min_h, user_max_h = usedAtoms.get('H')

        query_obj = self.sql_db.session.query(func.max(CarbonHydrogen.C).label("max_c"), 
                        func.min(CarbonHydrogen.C).label("min_c"),
                        func.max(CarbonHydrogen.H).label("max_h"),
                        func.min(CarbonHydrogen.H).label("min_h"),
                        )


        database = query_obj.first()
        if database.max_c == user_max_c and database.min_c == user_min_c and database.max_h == user_max_h and database.min_h == user_min_h:   
            #all data is already available at the database
            pass

        else:

            current_count = self.sql_db.session.query(CarbonHydrogen.C).count()
            
            databaseCarbonHydrogen = self.sql_db.session.query(CarbonHydrogen).all()
            
            userCarbon = set(range(user_min_c, user_max_c + 1))
            userHydrogen = set(range(user_min_h, user_max_h + 1))
            
            carbon_hydrogen_objs_database = {}
            for obj in databaseCarbonHydrogen:
                
                str_data = "C:{},H:{}".format(obj.C, obj.H)
                carbon_hydrogen_objs_database[str_data] = str_data

            carbon_hydrogen_objs_to_create = {'even': {}, 'odd': {}}
            
            list_ch_obj_to_add = list()
            i = 0
            for comb in itertools.product(userCarbon, userHydrogen):
                
                C  = comb[0]
                H =  comb[1]
                data = {"C":C,
                       "H":H,
                }

                data_insert = {"C":C,
                       "H":H,
                }

                str_data = "C:{},H:{}".format(C,H)
                
                if not str_data in carbon_hydrogen_objs_database.keys():
                    
                    label = 'even' if comb[1]%2 == 0 else 'odd'
                    data["mass"] = (C * Atoms.atomic_masses.get('C')) + (H * Atoms.atomic_masses.get('H'))
                    data["dbe"] = C - (H/2) + 1
                    data["id"] = i + current_count + 1
                    data_insert["id"] = i + current_count + 1
                    i = i + 1 
                    carbon_hydrogen_objs_to_create[label][str_data] = data
                    
                    list_ch_obj_to_add.append(data_insert)

            if list_ch_obj_to_add:
                # insert carbon hydrogen objs
                list_insert_chunks = chunks(list_ch_obj_to_add, self.sql_db.chunks_count)
                for insert_chunk in  list_insert_chunks:   
                    insert_query = CarbonHydrogen.__table__.insert().values(insert_chunk)
                    self.sql_db.session.execute(insert_query)
                self.sql_db.session.commit()    
            
                
                list_molecular_form= list()
                for classe_obj in existing_classes_objs:

                    classe_dict = classe_obj.to_dict()  
                    classe_mass = self.calc_mz(classe_dict)
                    classe_dbe = self.calc_dbe_class(classe_dict)

                    odd_even_label = self.get_h_odd_or_even(classe_dict)

                    ch_datalist = carbon_hydrogen_objs_to_create.get(odd_even_label).values()

                    for ch_dict in ch_datalist:
                        mass = ch_dict.get('mass') + classe_mass
                        dbe = ch_dict.get('dbe') + classe_dbe

                        if settings.min_mz <= mass <= settings.max_mz:
                
                            if settings.min_dbe <= dbe <= settings.max_dbe:
                                
                                list_molecular_form.append( {"heteroAtoms_id":classe_obj.id, 
                                        "carbonHydrogen_id":ch_dict.get('id'), 
                                        "mass":mass, "DBE":dbe})

                list_insert_chunks = chunks(list_molecular_form, self.sql_db.chunks_count)
                for insert_chunk in  list_insert_chunks:   
                    insert_query = MolecularFormulaLink.__table__.insert().values(insert_chunk)
                    self.sql_db.session.execute(insert_query)
                self.sql_db.session.commit()    
            
                
    @timeit
    def runworker(self, molecular_search_settings):
        
        classes_list, class_to_create, existing_classes_objs = self.check_database_get_class_list(molecular_search_settings)
        
        settings = MolecularLookupDictSettings()
        settings.usedAtoms = deepcopy(molecular_search_settings.usedAtoms)
        settings.url_database = molecular_search_settings.url_database
        settings.db_jobs = molecular_search_settings.db_jobs

        self.add_carbonsHydrogens(settings, existing_classes_objs)
        
        if class_to_create:
            
            settings = MolecularLookupDictSettings()
            settings.usedAtoms = deepcopy(molecular_search_settings.usedAtoms)
            settings.url_database = molecular_search_settings.url_database
            settings.db_jobs = molecular_search_settings.db_jobs
            
            self.sql_db.session.commit()
            odd_ch_obj = self.get_carbonsHydrogens(settings,'odd')
            self.odd_ch_id = [obj.id for obj in odd_ch_obj]
            self.odd_ch_dict = [{'C':obj.C, 'H':obj.H} for obj in odd_ch_obj]
            self.odd_ch_mass = [obj.mass for obj in odd_ch_obj]
            self.odd_ch_dbe = [obj.dbe for obj in odd_ch_obj]
            
            even_ch_obj = self.get_carbonsHydrogens(settings, 'even')
            self.even_ch_id = [obj.id for obj in even_ch_obj]
            self.even_ch_dict = [{'C':obj.C, 'H':obj.H} for obj in even_ch_obj]
            self.even_ch_mass = [obj.mass for obj in even_ch_obj]
            self.even_ch_dbe = [obj.dbe for obj in even_ch_obj]

            all_results= list()
            for class_tuple in tqdm(class_to_create):
                
                results = self.populate_combinations(class_tuple, settings)
                all_results.extend(results)
                if settings.db_jobs == 1: 
                    #if len(all_results) >= self.sql_db.chunks_count:
                        list_insert_chunks = list(chunks(results, self.sql_db.chunks_count))
                        for chunk in list_insert_chunks:
                            insert_query = MolecularFormulaLink.__table__.insert().values(chunk)
                            self.sql_db.session.execute(insert_query)
                        #all_results = list()
            self.sql_db.session.commit()
            # each chunk takes ~600Mb of memory, so if using 8 processes the total free memory needs to be 5GB
            if settings.db_jobs > 1: 
                list_insert_chunks = list(chunks(all_results, self.sql_db.chunks_count))
                print( "Started database insert using {} iterations for a total of {} rows".format(len(list_insert_chunks), len(all_results)))
                worker_args = [(chunk, settings.url_database) for chunk in list_insert_chunks]
                p = multiprocessing.Pool(settings.db_jobs)
                for class_list in tqdm(p.imap_unordered(insert_database_worker, worker_args)):
                    pass
                p.close()
                p.join()
        
        return classes_list
    
    def get_classes_in_order(self, molecular_search_settings):
        ''' structure is 
            ('HC', {'HC': 1})'''
        
        usedAtoms = deepcopy(molecular_search_settings.usedAtoms)
        
        usedAtoms.pop("C")
        usedAtoms.pop("H")

        min_n, max_n = usedAtoms.get('N') if usedAtoms.get('N') else (0,0)
        min_o, max_o = usedAtoms.get('O') if usedAtoms.get('O') else (0,0)
        min_s, max_s = usedAtoms.get('S') if usedAtoms.get('S') else (0,0)
        min_p, max_p = usedAtoms.get('P') if usedAtoms.get('P') else (0,0)

        possible_n = [n for n in range(min_n, max_n + 1)]
        possible_o = [o for o in range(min_o, max_o + 1)]
        possible_s = [s for s in range(min_s, max_s + 1)]
        possible_p = [p for p in range(min_p, max_p + 1)]
        
        atoms_in_order = ['N', 'O', 'S', 'P']

        classe_in_order = {}

        all_atoms_tuples = itertools.product(possible_n, possible_o,
                                            possible_s, possible_p)
        
        for atom in atoms_in_order:
            usedAtoms.pop(atom, None)
        
        for selected_atom, min_max_tuple in usedAtoms.items():
            
            min_x = min_max_tuple[0]
            max_x = min_max_tuple[1]
            

            possible_x = [x for x in range(min_x, max_x + 1)]

            all_atoms_tuples = itertools.product(all_atoms_tuples, possible_x)
            all_atoms_tuples = [all_atoms_combined[0] + (all_atoms_combined[1],) for all_atoms_combined in
                                all_atoms_tuples]
            atoms_in_order.append(selected_atom)
        
        for all_atoms_tuple in all_atoms_tuples:

            classe_str = ''
            classe_dict = {}
            
            for each_atoms_index, atom_number in enumerate(all_atoms_tuple):
                
                if atom_number != 0:
                    classe_dict[atoms_in_order[each_atoms_index]] = atom_number
            
            if not classe_dict:
                classe_in_order['HC'] = {"HC": ""}
                continue

            classe_str =json.dumps(classe_dict)
            
            if len(classe_str) > 0:
                
                classe_in_order[classe_str] =  classe_dict
        
        classe_in_order_dict = self.sort_classes(atoms_in_order, classe_in_order)
        
        return classe_in_order_dict

    @staticmethod
    def sort_classes( atoms_in_order, combination_dict) -> [str]: 
        #ensures atoms are always in the order defined at atoms_in_order list
        join_dict_classes = dict()
        atoms_in_order =  ['N','S','P','O'] + atoms_in_order[4:] + ['HC']
        
        sort_method = lambda atoms_keys: [atoms_in_order.index(atoms_keys)] 
        for class_str, class_dict in combination_dict.items():
            
            sorted_dict_keys = sorted(class_dict, key = sort_method)
            class_dict = { atom: class_dict[atom] for atom in sorted_dict_keys}
            class_str = json.dumps(class_dict)
            # using json for the new database, class 
            # class_str = ' '.join([atom + str(class_dict[atom]) for atom in sorted_dict_keys])
            join_dict_classes[class_str] =  class_dict
        
        return join_dict_classes

    @staticmethod
    def get_fixed_initial_number_of_hydrogen( min_h, odd_even):

        remaining_h = min_h % 2
        
        if odd_even == 'even':
            
            if remaining_h == 0: return remaining_h
            
            else: return remaining_h + 1    
        
        else:
            
            if remaining_h == 0: return remaining_h + 1
            
            else: return remaining_h    

    def calc_mz(self, datadict, class_mass=0):
        
        mass = class_mass
        
        for atom in datadict.keys():
            
            if atom != 'HC':    
                
                mass = mass + Atoms.atomic_masses[atom]  *  datadict.get(atom)
            
        return mass 
        
    def calc_dbe_class(self, datadict):
            
            init_dbe = 0
            for atom in datadict.keys():

                if atom == 'HC':
                    continue  
                
                n_atom = int(datadict.get(atom))
                
                clean_atom = ''.join([i for i in atom if not i.isdigit()]) 
                
                valencia = MSParameters.molecular_search.used_atom_valences.get(clean_atom)
                
                if type(valencia) is tuple:
                    valencia = valencia[0]
                if valencia > 0:
                    #print atom, valencia, n_atom, init_dbe
                    init_dbe = init_dbe + (n_atom * (valencia - 2))
                else:
                    continue
                
            return (0.5 * init_dbe)
            
    def populate_combinations(self, classe_tuple, settings):
        
        ion_charge =  0
        
        class_dict = classe_tuple[1]
        odd_or_even = self.get_h_odd_or_even(class_dict)
        
        return self.get_mol_formulas(odd_or_even, classe_tuple, settings)
        
    def get_or_add(self, SomeClass, kw):
            
        obj = self.sql_db.session.query(SomeClass).filter_by(**kw).first()
        if not obj:
            obj = SomeClass(**kw)
        return obj
    
    
    def get_mol_formulas(self, odd_even_tag, classe_tuple, settings):
        
        class_str = classe_tuple[0]
        class_dict = classe_tuple[1]
        classe_id = classe_tuple[2]
        
        results = list()
        
        if 'HC' in class_dict:
            del class_dict['HC']
            
        class_dbe = self.calc_dbe_class(class_dict)    
        class_mass = self.calc_mz(class_dict)
        
        carbonHydrogen_mass = self.odd_ch_mass if odd_even_tag == 'odd' else self.even_ch_mass 
        carbonHydrogen_dbe = self.odd_ch_dbe if odd_even_tag == 'odd' else self.even_ch_dbe 
        carbonHydrogen_id = self.odd_ch_id if odd_even_tag == 'odd' else self.even_ch_id 
        
        for index, carbonHydrogen_obj in enumerate(carbonHydrogen_id):
            
            mass = carbonHydrogen_mass[index] + class_mass
            dbe =  carbonHydrogen_dbe[index] + class_dbe
    
            if settings.min_mz <= mass <= settings.max_mz:
                
                if settings.min_dbe <= dbe <= settings.max_dbe:
                    
                    molecularFormula=  {"heteroAtoms_id":classe_id, 
                            "carbonHydrogen_id":carbonHydrogen_id[index], 
                            "mass":mass, "DBE":dbe}
                    
                    results.append(molecularFormula)
        
        return results
        
        
    def get_h_odd_or_even(self, class_dict):

        
        HAS_NITROGEN = 'N' in class_dict.keys()
        HAS_PHOSPHORUS = 'P' in class_dict.keys()

        number_of_halogen = self.get_total_halogen_atoms(class_dict)

        if number_of_halogen > 0:

            TEM_HALOGEN = True

        else:

            TEM_HALOGEN = False

        if TEM_HALOGEN:

            remaining_halogen = number_of_halogen % 2

        else:

            remaining_halogen = 0

        if HAS_NITROGEN and HAS_PHOSPHORUS:

            number_of_n = class_dict.get('N') + class_dict.get('P')
            remaining_n = number_of_n % 2

        elif HAS_NITROGEN and not HAS_PHOSPHORUS:

            number_of_n = class_dict.get('N')
            remaining_n = number_of_n % 2

        elif HAS_PHOSPHORUS and not HAS_NITROGEN:

            number_of_n = class_dict.get('P')
            remaining_n = number_of_n % 2

        else:

            remaining_n = -1

        if remaining_n > 0.0:
            if HAS_NITROGEN or HAS_PHOSPHORUS:

                if TEM_HALOGEN:
                    if remaining_halogen == 0:
                        return 'odd'
                    else:
                        return 'even'
                else:
                    return 'odd'

        elif remaining_n == 0.0:

            if HAS_NITROGEN or HAS_PHOSPHORUS:

                if TEM_HALOGEN:
                    if remaining_halogen == 0:
                        return 'even'
                    else:
                        return 'odd'
                else:
                    return 'even'

        else:

            if TEM_HALOGEN:
                if remaining_halogen == 0:
                    return 'even'
                else:
                    return 'odd'
            else:
                return 'even'

    @staticmethod
    def get_total_halogen_atoms(class_dict):

            atoms = ['F', 'Cl', 'Br', 'I']

            total_number = 0
            
            for atom in atoms:

                if atom in class_dict.keys():

                    total_number = total_number + class_dict.get(atom)
            
            return total_number    
