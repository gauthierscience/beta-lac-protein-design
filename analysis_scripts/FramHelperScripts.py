############################
#                          #
#                          #
#     COMMON FUNCTIONS     #
#           FOR            #
#        FRAM ET AL        #
#                          #
#                          #
############################

import pandas as pd
import numpy as np

class FramHelperFunctions:
	'''
	A set of helper functions for analysis of the data
	from Fram et al., 2023
	'''

	def __init__(self, data_dir):
		self.data_dir = data_dir
		self.data_filename = '{0}/{1}'.format(data_dir, 'Raw Data.xlsx')
		self.model_filename = '{0}/{1}'.format(data_dir, 'BLAT_ECOLX_1_b0.5.model')
		self.msa_filename = '{0}/{1}'.format(data_dir, 'BLAT_ECOLX_1_b0.5.a2m')
		self.pdb_1xpb_filename = '{0}/{1}'.format(data_dir, '1xpb.pdb')

		self.sample_names_df = pd.read_excel(
		    self.data_filename, sheet_name='sample_order',
		)[['sample_name', 'synonyms']].astype(str).assign(
			synonyms = lambda df: df.synonyms.apply(
				lambda synonyms: synonyms.split(',')
            )
        )
	
	#
	#
	# FILENAMES
	#
	#
	def get_data_filename(self):
		return self.data_filename


	#
	#
	# SEQUENCES
	#
	#
	def get_sequences_df(self, only_tested=True):
		if not only_tested:
			return self._parseSequences(only_tested)

		if hasattr(self, 'sequences_df'): return self.sequences_df
		self.sequences_df = self._parseSequences()
		return self.sequences_df
	
	def _parseSequences(self, only_tested=True):
		sheetname = 'tested_seqs'
		if not only_tested: sheetname = 'all_generated_design_seqs'
        
		toreturn = self.add_manuscript_name_to_df(
			pd.read_excel(
				self.data_filename, sheet_name=sheetname
			), 
			synonym_column='sample_id', 
			new_column='manuscript_name'
		).drop('sample_id', axis=1)

		#load mutations
		wt_full_seq = toreturn[
			toreturn.manuscript_name=='WT TEM-1'
		].iloc[0].full_sequence

		numbering_df = self.get_numbering_df()
		def parseMutations(full_seq):
			muts = []
			for idx, seq_aa in enumerate(full_seq):
				if wt_full_seq[idx] != seq_aa:
					muts.append(
						'{0}{1}{2}'.format(
							wt_full_seq[idx], 
							int(
								numbering_df[
									numbering_df.model_full_seq_idx==idx
								].iloc[0].model_num
							),
							seq_aa
						)
					)
			return muts

		toreturn['model_mutations'] = toreturn.full_sequence.apply(
			lambda seq: parseMutations(seq)
		)

		#move manuscript name to the front of the dataframe
		col = toreturn.pop('manuscript_name')
		toreturn.insert(0, col.name, col)
		return toreturn


	def get_natural_msa_df(self, remove_unaligned_positions=True):
		'''
		parse in the natural multiple sequence alignment used to generate the model
	       
	    *requires !pip install biopython
		'''
		if hasattr(self, 'msa_df'): return self.msa_df
		from Bio import SeqIO

		with open(self.msa_filename) as fasta_file:  # Will close handle cleanly
			toreturn = {
				'seq_id': [],
				'seq': []
			}
			numbering_df = self.get_numbering_df()
			unaligned_positions = list(
				numbering_df[
					numbering_df.model_num.isna()
				].model_full_seq_idx
			)

			for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
				toreturn['seq_id'].append(seq_record.id)

				seq = ''.join( seq_record.seq )
				if remove_unaligned_positions:
					seq = ''.join([
						char for idx, char in enumerate(seq) if idx not in unaligned_positions
					])

				toreturn['seq'].append(seq)
			
			return pd.DataFrame(toreturn)

	#
	#
	# SAMPLE / MANUSCRIPT NAMING
	#
	#
	def get_sample_names_df(self):
		'''
		get the raw sample names dataframe
		'''
		return self.sample_names_df

	def add_manuscript_name_to_df(self, df, synonym_column, new_column='manuscript_name'):
		'''
		add a column to an arbritrary dataframe that contains the manuscript name for each
		row, translating from a synonym column
		raises an error if any of the synonyms are not found
		'''
		df = df.copy()

		selection_mask = lambda synonym: self.sample_names_df[(
			(self.sample_names_df.sample_name == str(synonym)) | 
			(self.sample_names_df.synonyms.apply(lambda x: str(synonym) in x))
		)]

		try:
			df[new_column] = df[synonym_column].apply(
				lambda synonym: selection_mask(synonym).iloc[0].sample_name
			)
		except IndexError:
			for synonym in df[synonym_column]:
				if len(selection_mask(synonym)) < 1:
					raise Exception('Invalid synonym found: {0}'.format(synonym))
			raise Exception('Unknown Error Occurred!')
		return df
	
	def get_sample_order(self, valid_manuscript_names = None):
		'''
		get an ordered list of samples as it appears in the manuscript, dropping
		ay samples not included in 'valid_manuscript_names' (if provided)
		'''
		toreturn = []
		for manuscript_name in self.sample_names_df.sample_name:
			if valid_manuscript_names is None or manuscript_name in list(valid_manuscript_names):
				toreturn.append(manuscript_name)
		return toreturn

	#
	#
	# GENERAL NUMBERING HELPER FUNCTIONS - PDB, MODEL, DMS
	#
	#
	def get_numbering_df(self):
		'''
		returns a dataframe that contains the WT TEM-1 sequence aligned
		to the model, dms, and pdb 1xpb - useful for deciphering the
		numbering differences
		'''
		if hasattr(self, 'numbering_df'): return self.numbering_df
		self.numbering_df = pd.read_excel(
			self.data_filename, sheet_name='numbering',
			header=13
		)
		self.numbering_df.model_num = self.numbering_df.model_num.astype(pd.Int64Dtype())
		return self.numbering_df
	
	def add_pdb_num_from_model_num(self, df, model_num_col, new_pdb_col):
		'''
		add a column to the dataframe that contains the pdb numbering from
		the model numbering
		'''
		numbering_df = self.get_numbering_df()
		toreturn = df.copy()
		toreturn[new_pdb_col] = toreturn[model_num_col].apply(
			lambda model_num: numbering_df[numbering_df.model_num == model_num].iloc[0].pdb_num
		)
		return toreturn
	

	def add_model_num_from_pdb_num(self, df, pdb_num_col, new_model_col):
		'''
		add a column to the dataframe that contains the model numbering from
		the pdb numbering
		'''
		numbering_df = self.get_numbering_df()
		toreturn = df.copy()
		toreturn[new_model_col] = toreturn[pdb_num_col].apply(
			lambda pdb_num: numbering_df[numbering_df.pdb_num == pdb_num].iloc[0].model_num
		).astype(pd.Int64Dtype())
		return toreturn
	
	def add_model_num_from_dms_num(self, df, dms_num_col, new_model_num_col):
		'''
		add a column to the dataframe that contains the model numbering from
		the dms numbering
		'''
		numbering_df = self.get_numbering_df()
		toreturn = df.copy()
		toreturn[new_model_num_col] = toreturn[dms_num_col].apply(
			lambda dms_num: numbering_df[numbering_df.dms_num == dms_num].iloc[0].model_num
		).astype(pd.Int64Dtype())
		return toreturn

	#
	#
	# MODEL
	#
	# REQUIRES: !pip install evcouplings
	#
	#
	def get_model(self):
		'''
		get the model
		'''
		from evcouplings.couplings import CouplingsModel
		if hasattr(self, 'model'): return self.model
		self.model = CouplingsModel(self.model_filename)
		return self.model

	def get_model_mutation_effect_table(self):
		'''
		get a dataframe of all singles for both the full epistatic model
		prediction as well as the independent model prediction
		'''
		model = self.get_model()

		from evcouplings.mutate.calculations import single_mutant_matrix, predict_mutation_table

		#mutational effect prediction
		df_evh_singles = single_mutant_matrix(
			model, output_column='effect_prediction_epistatic'
		)

		#add independent model predictions
		df_evh_singles = predict_mutation_table(
			model.to_independent_model(), 
			df_evh_singles,
			output_column = 'effect_prediction_independent'
		)
		return df_evh_singles
	
	#
	#
	# PDB of WT TEM-1 - 1XPB
	# https://www.rcsb.org/3d-view/1XPB
	#
	# requires !pip install biopython
	#
	#
	def renumber_1xpb_pdb_chain(self, pdb_chain):
		'''
		the 1xpb pdb file has two position 51s: 51 and 51a and is missing
		position 57 - this function returns a corrected pdb chain where
		positions 51A->56 become 52-57. All other residues are unaltered.

		resi.id: the biopython pdb parser value available for each 
		         residue, acquired by calling "resi.id"
		'''
		for idx in reversed(range(52, 57)):
			pdb_chain[idx].id = (' ', idx+1, ' ')
		pdb_chain[(' ',51,'A')].id = (' ', 52, ' ')
		return pdb_chain


	def get_1xpb_pdb_chain_polypeptide(self):
		'''
		load 1xpb and fix numbering 
		'''
		from Bio.PDB import PDBParser

		pdb_parser = PDBParser(QUIET=True)
		pdb_structure = pdb_parser.get_structure('1xpb', self.pdb_1xpb_filename)
		pdb_model = pdb_structure[0] 
		pdb_chain = pdb_model['A']

		pdb_chain = self.renumber_1xpb_pdb_chain(pdb_chain)
		return pdb_chain
	
	def get_1xpb_interactions(self, angstrom_cutoff=5):
		'''
		This function returns a dataframe of residue by residue interactions in
		the 1xpb file. Interaction is defined as having any atom in residue 1 
		within a distance of 5 angstroms of any any atom in residue 2
		'''
		from Bio.PDB.Polypeptide import protein_letters_3to1
		
		pdbchain = self.get_1xpb_pdb_chain_polypeptide()
		pdb_interactions = {
			'pdb_i': [], 'pdb_j': [], 
			'pdb_i_aa': [], 'pdb_j_aa': [], 
			'min_atom_distance': []
		}

		all_aa_residues = [res for res in pdbchain if 'CA' in res]
		for r1_idx, r1 in enumerate(all_aa_residues):
			r1_num = r1.id[1]
			r1_aa = protein_letters_3to1[r1.get_resname()]
			for r2_idx, r2 in enumerate(all_aa_residues[r1_idx:]):
				r2_num = r2.id[1]
				r2_aa = protein_letters_3to1[r2.get_resname()]
				if r1 != r2:
					distance = None

					min_distance = 1000000000
					for atom1 in r1:
						for atom2 in r2:
							distance = np.linalg.norm(atom1.coord - atom2.coord)
							if distance < min_distance: min_distance = distance
					distance = min_distance
					if distance <= angstrom_cutoff:
						pdb_interactions['pdb_i'].append(r1_num)
						pdb_interactions['pdb_j'].append(r2_num)
						pdb_interactions['pdb_i_aa'].append(r1_aa)
						pdb_interactions['pdb_j_aa'].append(r2_aa)
						pdb_interactions['min_atom_distance'].append(min_distance)
		
		toreturn = pd.DataFrame(pdb_interactions)
		toreturn = self.add_model_num_from_pdb_num(
			toreturn, 'pdb_i', 'model_i'
		)
		toreturn = self.add_model_num_from_pdb_num(
			toreturn, 'pdb_j', 'model_j'
		)
		return toreturn

	#
	#
	# DSSP of WT TEM-1 - 1XPB
	# downloaded from MRS - https://mrs.cmbi.umcn.nl/
	#
	# requires !pip install biopython
	#
	def get_dssp(self):
		if hasattr(self, 'dssp'): return self.dssp.copy()
		self.dssp = self._parseDssp()
		return self.dssp

	def _parseDssp(self):
		'''
        This function takes the complex output of the BioPython DSSP parser
        and converts it into a pandas dataframe.
        '''
		from Bio.PDB.DSSP import make_dssp_dict
		
		#surface_area from http://prowl.rockefeller.edu/aainfo/volume.htm
		residue_surface_area = {
			'A': 115, 'R': 225,'D': 150, 'N': 160, 'C': 135, 'E': 190, 'Q': 180, 
			'G': 75, 'H': 195, 'I': 175, 'L': 170, 'K': 200, 'M': 185, 'F': 210,
			'P': 145, 'S': 115, 'T': 140, 'W': 255, 'Y': 230, 'V': 155,
		}

		toreturn = {
            'position': [],'chain': [],
            'aa': [], 'ss': [], 'acc': [], 'phi': [], 'psi': [],
            'N-H-->O [1]': [], 'N-H-->O [2]': [],
            'O-->N-H [1]': [], 'O-->N-H [2]': [],
            'relative_acc': [], 'aa_surface_area': [],
        }

		raw_dssp_dict = make_dssp_dict(self.data_dir+'/1xpb_numbering_corrected.dssp')
		for residue_key in raw_dssp_dict[0].keys(): #residue_key: ('chainid', ('?', position, '?'))
			data = raw_dssp_dict[0][residue_key]
			toreturn['position'].append(
				residue_key[1][1]
			)
			toreturn['chain'].append(residue_key[0])

			toreturn['aa'].append(data[0])
			toreturn['ss'].append(data[1])
			toreturn['acc'].append(data[2])
			toreturn['phi'].append(data[3])
			toreturn['psi'].append(data[4])

			toreturn['N-H-->O [1]'].append([data[6], data[7]])
			toreturn['N-H-->O [2]'].append([data[10], data[11]])
			toreturn['O-->N-H [1]'].append([data[8], data[9]])
			toreturn['O-->N-H [2]'].append([data[12], data[13]])

			#calculate relative solvent accessable area
			AA = data[0]
			if AA == 'a': AA = 'C' #the lower case a represents cystine
			toreturn['aa_surface_area'].append( residue_surface_area[AA] )
			toreturn['relative_acc'].append(
				data[2] / residue_surface_area[AA]
			)
		return pd.DataFrame(toreturn)

	#
	#
	# STIFFLER ET AL DMS
	# doi: 10.1016/j.cell.2015.01.035. 
	# PMID: 25723163.
	#
	#
	def get_dms_df(self):
		if hasattr(self, 'dms_df'): return self.dms_df.copy()
		self.dms_df = self._parseStiffler()
		return self.dms_df
	
	def _parseStiffler(self):
		toreturn = pd.merge(
			pd.read_excel(
				self.data_filename, sheet_name='stiffler_dms_rep1',
			).melt(
				id_vars=['Mutation'], var_name='amp_conc'
			),
			pd.read_excel(
				self.data_filename, sheet_name='stiffler_dms_rep2',
			).melt(
				id_vars=['Mutation'], var_name='amp_conc'
			),
			on=['Mutation', 'amp_conc']
		).rename(
			columns={
				'Mutation': 'stiffler_mutation',
				'value_x': 'rep1_value', 
				'value_y': 'rep2_value'
			}
		)

		toreturn['mut_from'] = toreturn.stiffler_mutation.apply(lambda mut: mut[0])
		toreturn['stiffler_position'] = toreturn.stiffler_mutation.apply(lambda mut: int(mut[1:-1]))
		toreturn['mut_to'] = toreturn.stiffler_mutation.apply(lambda mut: mut[-1])

		toreturn = self.add_model_num_from_dms_num(toreturn, 'stiffler_position', 'model_pos')
		toreturn['model_mutation'] = toreturn.apply(lambda row:
			'{0}{1}{2}'.format(row['mut_from'], row['model_pos'], row['mut_to']), axis=1
		)
		toreturn = toreturn.dropna(subset='model_pos') #remove positions not aligned in the model

		#
		# add design names to rows in which designs have mutations
		#
		def get_mut_seqlist_df(seqs, manuscript_names_col):
			'''
			returns a dataframe with columns:
				'model_mutation' - a string
				manuscript_names_col - from the parameter - an array
			'''
			return seqs.explode('model_mutations')[
				['manuscript_name', 'is_design', 'functional', 'model_mutations']
			].dropna().groupby(
				'model_mutations'
			).agg(
				{'manuscript_name': lambda x: x.tolist()}
			).reset_index().rename(
				columns={
					'model_mutations': 'model_mutation',
					'manuscript_name': manuscript_names_col,
				}
			)

		seqs_df = self.get_sequences_df()
		toreturn = pd.merge(
			toreturn, 
			get_mut_seqlist_df(
				seqs_df[seqs_df.is_design == 0], 'ctrls_with_mut'
			),
			how='left', on='model_mutation', 
		).assign(
			ctrls_with_mut = lambda df: [
				[] if x is np.NaN else x for x in df.ctrls_with_mut
			]
		)

		toreturn = pd.merge(
			toreturn, 
			get_mut_seqlist_df(
				seqs_df[((seqs_df.is_design == 1) & (seqs_df.functional == 1))], 
				'func_designs_with_mut'
			),
			how='left', on='model_mutation', 
		).assign(
			func_designs_with_mut = lambda df: [
				[] if x is np.NaN else x for x in df.func_designs_with_mut
			]
		)

		toreturn = pd.merge(
			toreturn, 
			get_mut_seqlist_df(
				seqs_df[((seqs_df.is_design == 1) & (seqs_df.functional == 0))], 
				'nonfunc_designs_with_mut'
			),
			how='left', on='model_mutation', 
		).assign(
			nonfunc_designs_with_mut = lambda df: [
				[] if x is np.NaN else x for x in df.nonfunc_designs_with_mut
			]
		)

		return toreturn