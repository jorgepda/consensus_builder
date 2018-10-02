#!/data/perezde/dna_barcoding/dna_barcoding_env3/bin/python

import os
import inspect
import sys
import re
import getopt
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

accepted_formats = ["fasta", "fastq"]
default_format = "fasta"

# Script that finds ab1 pairs, trims them, builds a consensus, and retrims the consensus.
# All ab1 files must be contained in the same directory.
## Usage:
# find_pairs("/path/to/directory)


def longest_common_substring(str1, str2):
	m = [[0] * (1 + len(str2)) for i in range(1 + len(str1))]
	longest, x_longest = 0, 0
	for x in range(1, 1 + len(str1)):
		for y in range(1, 1 + len(str2)):
			if str1[x - 1] == str2[y - 1]:
				m[x][y] = m[x - 1][y - 1] + 1
				if m[x][y] > longest:
					longest = m[x][y]
					x_longest = x
			else:
				m[x][y] = 0
	return str1[x_longest - longest: x_longest]	

def trim(seq, cutoff=0.05, window=20, save=True, output=None, file_format=default_format, filter_Ns=False):
	start = False
	trim_start = 0
	
	if save and (output is None or not os.path.isdir(output)):
		print("Please provide a valid directory. Exiting program...")
		return None  	

	score_list = [cutoff - (10 ** (qual / -10.00)) for qual in seq.letter_annotations['phred_quality']]
	cummul_score = [score_list[0]]

	for i in range(1, len(score_list)):
		score = cummul_score[i - 1] + score_list[i]
		
		if score < 0:
			cummul_score.append(0)
		else:
			cummul_score.append(score)
			if not start:
				trim_start = i
				start = True
		
		if filter_Ns and seq.seq[i] == 'N' and start and (i - trim_start < window):
			start = False		
	
	trim_finish = cummul_score.index(max(cummul_score))

	end_window = 0
	
	while end_window < window:
		if trim_finish - end_window > 0 and seq.seq[trim_finish - end_window] == 'N':
			trim_finish = trim_finish - end_window - 1
			end_window = 0			

		end_window += 1
	
	trimmed_sequence = seq[trim_start:trim_finish]

	if save and output is not None and os.path.isdir(output):
		if file_format not in accepted_formats:
			print(file_format + " is not an accepted format. Will save as fasta.")
			file_format = default_format

		if output[len(output) - 1] != '/':
			output += '/'
		
		if trimmed_sequence.name != '':
			file_name = trimmed_sequence.name + "_trimmed"
		else:
			file_name = "unnamed_trimmed_sequence"

		file_name = output + file_name + "." + file_format		

		try:
			SeqIO.write(trimmed_sequence, file_name, file_format)
		except:
			e = sys.exc_info()[1]
			print("Trimmed file " + file_name + " could no be created with error: " + e)		
	
	return trimmed_sequence

def align(seq_1, seq_2, gap_open = 50, gap_extend = 5):
	aligner = Align.PairwiseAligner()

	aligner.open_gap_score = -1 * gap_open
	aligner.extend_gap_score = -1 * gap_extend
#	aligner.target_end_gap_score = 0
#	aligner.query_end_gap_score = 0
	
	if len(seq_1.seq) > 0 and len(seq_2.seq) > 0:
		alignment = aligner.align(seq_1.seq, seq_2.seq.reverse_complement())
	elif len(seq_1.seq) > 0:
		alignment = [seq_1.seq]
	elif len(seq_2.seq) > 0:
		alignment = [seq_2.seq.reverse_complement()]
	else:
		alignment = [None]

	return alignment[0]	

def get_consensus(alignment, forward_record, reverse_record, save=True, output=None, file_format=default_format):
	consensus = []
	phred_scores = []
	forward_pointer = -1
	reverse_pointer = -1
	forward_phred_scores = forward_record.letter_annotations['phred_quality']
	reverse_phred_scores = reverse_record.reverse_complement().letter_annotations['phred_quality']

	if save and (output is None or not os.path.isdir(output)):
		print("Please provide a valid directory. Exiting program...")
		return None  	

	forward = str(alignment).splitlines()[0]
	path = str(alignment).splitlines()[1]
	reverse = str(alignment).splitlines()[2]

	for i in range(0, len(path)):
		if path[i] == '|':
			consensus.append(forward[i])
			forward_pointer += 1
			reverse_pointer += 1
			phred_scores.append(forward_phred_scores[forward_pointer] if forward_phred_scores[forward_pointer] > reverse_phred_scores[reverse_pointer] else reverse_phred_scores[reverse_pointer])
		elif path[i] == '-':
			if forward[i] != '-':
				consensus.append(forward[i])
				forward_pointer += 1
				phred_scores.append(forward_phred_scores[forward_pointer])
			else:
				consensus.append(reverse[i])
				reverse_pointer += 1
				phred_scores.append(reverse_phred_scores[reverse_pointer])
		else: # there's a conflict
				forward_pointer += 1
				reverse_pointer += 1
				if forward_phred_scores[forward_pointer] > reverse_phred_scores[reverse_pointer]:
					consensus.append(forward[i])
					phred_scores.append(forward_phred_scores[forward_pointer])
				else:
					consensus.append(reverse[i])
					phred_scores.append(reverse_phred_scores[reverse_pointer])
	
	consensus = ''.join(consensus) 
	new_id = longest_common_substring(forward_record.name, reverse_record.name)
	consensus_record = SeqRecord(Seq(consensus, IUPAC.unambiguous_dna), id=new_id, name=new_id)
	consensus_record.letter_annotations['phred_quality'] = phred_scores
	description = "Taxonomy unknown. Consensus for reads " + forward_record.name + " " + reverse_record.name
	consensus_record.description = description

	if save and output is not None and os.path.isdir(output):
		if file_format not in accepted_formats:
			print(file_format + " is not an accepted format. Will save as fasta.")
			file_format = default_format

		if output[len(output) - 1] != '/':
			output += '/'
		
		file_name = output + new_id + "." + file_format 

		try:
			SeqIO.write(consensus_record, file_name, file_format)
		except:
			e = sys.exc_info()[1]
			print("Consensus file " + file_name + " could no be created with error: " + e)		


	return consensus_record

def read_abi(file_path):
	if file_path is not None and os.path.isfile(file_path):
		read = SeqIO.read(file_path, 'abi')
		return read
	else:
		print("Please enter a valid file path.")
		return None

def build_consensus(forward=None, reverse=None, save=True, output=None, file_format="fasta", trimming_window=20):
	
	if file_format not in accepted_formats:
		print(file_format + " is not an accepted format. Will save as fasta.")
		file_format = default_format

	if save and (output is None or not os.path.isdir(output)):
		print("Please provide a valid directory. Exiting program...")
		return None  	
	if forward is None and reverse is None:
		print("Please enter at least one file. Exiting program...")
		return None
	elif forward is not None and reverse is not None:
		forward = read_abi(forward)
		reverse = read_abi(reverse)
		forward_trimmed = trim(forward, save=False)
		reverse_trimmed = trim(reverse, save=False)
		if len(forward_trimmed.seq) > 0 and len(reverse_trimmed.seq) > 0:
			alignment = align(forward_trimmed, reverse_trimmed)
			consensus_sequence = get_consensus(alignment, forward_trimmed, reverse_trimmed, save=False)
			consensus_trimmed = trim(consensus_sequence, cutoff=0.01, save=save, window=trimming_window, output=output, file_format=file_format, filter_Ns = True)
		elif len(forward_trimmed.seq) > 0:
			consensus_trimmed = trim(forward_trimmed, cutoff=0.01, save=save, window=trimming_window, output=output, file_format=file_format, filter_Ns = True)
		elif len(reverse_trimmed.seq) > 0:
			reverse_trimmed.seq = reverse_trimmed.seq.reverse_complement()
			reverse_trimmed.letter_annotations['phred_quality'].reverse()
			consensus_trimmed = trim(reverse_trimmed, cutoff=0.01, save=save, window=trimming_window, output=output, file_format=file_format, filter_Ns = True)
		else:
			print("Files " + forward.name + " " + reverse.name + " did not pass quality parameters. No consensus file has been made.")
			return None

	elif forward is not None:
		forward = read_abi(forward)
		consensus_trimmed = trim(forward, cutoff=0.01, save=save, window=trimming_window, output=output, file_format=file_format, filter_Ns=True)
	else:
		# trim reverse only
		reverse = read_abi(reverse)
		reverse.seq = reverse.seq.reverse_complement()
		reverse.letter_annotations['phred_quality'].reverse_complement()
		consensus_trimmed = trim(reverse, cutoff=0.01, save=save, window=trimming_window, output=output, file_format=file_format, filter_Ns=True)

	if not save:
		print(consensus_trimmed.format(file_format))		

def quality_selection(seq, directory=None, trimming_window=10):
	if seq is None:
		print("Function takes array of sequences. None received. Exiting program...")
		return None
	
	if directory is not None and directory[len(directory) - 1] != '/':
		directory += '/'

	if len(seq) == 0:
		return None
	if len(seq) == 1:
		max_read = seq[0]
		if directory is not None:
			max_read = directory + max_read
	else:
		read = read_abi(directory + seq[0])
		if len(read.seq) > 2*trimming_window:
			current_max = sum(read.letter_annotations['phred_quality'][trimming_window:len(read.seq)-trimming_window])
		else:
			current_max = 0
		max_read = seq[0]

		for i in range(1, len(seq)):
			read = read_abi(directory + seq[i])
			if len(read.seq) > 2*trimming_window:
				current_qual = sum(read.letter_annotations['phred_quality'][trimming_window:len(read.seq)-trimming_window])
			else:
				current_qual = 0

			if current_qual >= current_max:
				current_max = current_qual
				max_read = seq[i]

		if directory is not None:
			max_read = directory + max_read

	return max_read
			

def find_pairs(directory, latest_resequence = True, save=True, output=None, file_format=default_format, trimming_window=20):
	if directory is None or not os.path.isdir(directory):
		print("Please provide a directory. Exiting execution...")
		return None

	if save and (output is None or not os.path.isdir(output)):
		print("Please provide a valid directory. Exiting program...")
		return None  	

	if directory[len(directory) - 1] != '/':
		directory += '/'
	
	abi_files = [file for file in os.listdir(directory) if re.search(r'[A-Z]{2,3}[-_]\d{2,3}.*?[FR].*?\.ab1', file)]
	abi_files.sort()
	searching_resequence = True
	forward_sequences = []
	reverse_sequences = []

	while len(abi_files)>0:
		current_file = abi_files.pop(0)
		match = re.search(r"^(.*?)F_", current_file)
		if match is not None:
			file_base = match.groups()[0]
			if len(abi_files) > 0:
				reverse_file = abi_files[0]
				match = re.search(r"^(.*?)R_", reverse_file)
			else:
				match = None
			if match is None: 
				forward_sequences.append(current_file)
			else:
				forward_sequences.append(current_file)
				if file_base == match.groups()[0]:
					reverse_file = abi_files.pop(0)
					reverse_sequences.append(reverse_file)
					while searching_resequence:
						if len(abi_files) > 0 and re.search(r"^" + file_base + r"(.*?)R_R_", abi_files[0]):
							reverse_file = abi_files.pop(0)
							reverse_sequences.append(reverse_file) 
						else:
							searching_resequence = False
				else:
					if re.search(r"^" + file_base, reverse_file):
						continue
		else:
			match = re.search(r"^(.*?)R_", current_file)
			if match is not None:
				file_base = match.groups()[0]
				reverse_sequences.append(current_file)
				while searching_resequence:
					if len(abi_files) > 0 and re.search(r"^" + file_base + r"R_R_", abi_files[0]):
						current_file = abi_files.pop(0)	
						reverse_sequences.append(current_file)
					else:
						searching_resequence = False
			else:
				print("File name " + current_file + " is not in the accepted format.")
				continue
	
		searching_resequence = True
		# print(forward_sequences)
		# print(reverse_sequences)

		if latest_resequence:
			if len(forward_sequences) > 0 and len(reverse_sequences) > 0:
				forward = directory + forward_sequences.pop()
				reverse = directory + reverse_sequences.pop()
				print("Files " + forward + " and " + reverse + " are paired.")
				consensus_sequence = build_consensus(forward=forward, reverse=reverse, save=save, output=output, file_format=file_format, trimming_window=trimming_window)
			elif len(forward_sequences) > 0:
				forward = directory + forward_sequences.pop()
				print("File " + forward + " is single ended.")
				build_consensus(forward=forward, reverse=None, save=save, output=output, file_format=file_format, trimming_window=trimming_window)
			else:
				reverse = directory + reverse_sequences.pop()
				print("File " + reverse + " is single ended.")
				build_consensus(forward=None, reverse=reverse, save=save, output=output, file_format=file_format, trimming_window=trimming_window)
		else:
			forward = quality_selection(forward_sequences, directory=directory)
			reverse = quality_selection(reverse_sequences, directory=directory)
			build_consensus(forward, reverse, save=save, output=output, file_format=file_format, trimming_window=trimming_window)

		print()
		forward_sequences = []
		reverse_sequences = []

def main(argv):
	# latest resequence
	# directory
	# save
	# target directory
	# format
	# trimming window
	input_path = None # d
	save = False # s
	latest_resequence = False # l
	output = None # t
	file_format = default_format #f
	trimming_window = 20 #w

	try:
		opts, args = getopt.getopt(argv,"hi:lo:f:w:",["help", "input=","latest-resequence","output=","file-format=","window="])
	except getopt.GetoptError:
		print("consensus.py -i <path> -o <path>")
		print("")
		sys.exit(2)
	for opt,arg in opts:
		if opt in ('-h', "--help"):
			print("consensus.py -i <path> -o <path>")
			print("")
			print("consensus.py:    script that given an input path where ab1 files are located, automatically pairs the files, aligns them, and creates a trimmed consensus. If no output path is given, consensus will be printed.")
			print("")
			print("Arguments:")
			print("-i, --input:\t\t\tpath where ab1 files are located")
			print("-l, --latest-resequence:\tselect newest resequence when aligning. If not, select based on quality. Default: False")
			print("-o, --output:\t\t\tpath where output files will be saved. If not provided, consensus will be printed. Default: None")
			print("-f, --file-format:\t\tfile format in which output wil be saved. Default: fasta. ")
			print("-w, --window:\t\t\tlength of window at read ends where Ns are not permitted. Default: 20")
			print("-h, --help:\t\t\tprints this screen and exits")
			sys.exit()
		elif opt in ("-i", "--input"):
			input_path = arg
		elif opt in ("-l", "--latest-resequence"):
			latest_resequence = True
		elif opt in ("-o", "--output"):
			output = arg
		elif opt in ("-f", "--file-format"):
			file_format = arg
		elif opt in ("-w", "--window"):
			trimming_window = int(arg)

	if input_path is None or not os.path.isdir(input_path):
		print("Please provide a valid path with ab1 files.")
		sys.exit()

	if output is not None:
		save = True 
		if not os.path.isdir(output):
			print("Please provide a valid path where consensus files will be saved.")
			sys.exit()

	if file_format not in accepted_formats:
		print(file_format + " not recognized. Using fasta instead.")
		file_format = default_format

	find_pairs(input_path, latest_resequence=latest_resequence, save=save, output=output, file_format=file_format, trimming_window=trimming_window)

if __name__ == "__main__":
	main(sys.argv[1:])
