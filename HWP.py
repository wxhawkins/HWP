import math
from functools import partial
import numpy as np
import re
import argparse
import time

from bokeh.models import HoverTool, Button, Label, LabelSet
from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh import palettes

from Bio.SubsMat import MatrixInfo
from Bio import SeqIO

#Establish argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-s", action = "store", type = str, dest = "seq_path", required = True)
parser.add_argument("-a", action = "store", type = str, dest = "align_path", default = None)
parser.add_argument("-r", action = "store", type = str, dest = "seq_range", default = None)
parser.add_argument("-c", action = "store_true", dest = "coiled_coil", default = False)
parser.add_argument("-f", action = "store_true", dest = "find", default = False)
args = parser.parse_args()


#Get primary sequence to be analyzed
try:
	seq = SeqIO.read(args.seq_path, "fasta").seq
	if len(seq) == 0:
		raise Exception
except:
	raise ValueError("Invalid protein sequence provided for analyzation.")

#Get and clean sequence alignment
align_seqs = args.align_path
if align_seqs is not None:
	align_seqs = []

	for fasta in SeqIO.parse(args.align_path, "fasta"):
		align_seqs.append(fasta.seq)

	align_seqs = np.array(align_seqs)

	gap_found = True
	while gap_found:
		gap_found = False
		for pos in range(len(align_seqs[0])):
			if align_seqs[0, pos] == "-":
				align_seqs = np.delete(align_seqs, pos, axis = 1)
				gap_found = True
				break


class Rule_set:
	def __init__(self, coiled_coil = False):
		#Establish alphabet
		self.all_aas = "WIFLCMVYPATHGSQNEDKR"

		#Set colors for each amino acid based on hydrophobicity (two aa's per color)
		pal = palettes.all_palettes['RdYlBu'][10]
		self.color_dict = {}
		for i, aa in enumerate(self.all_aas):
			color_ = pal[math.floor(i / 2)]
			self.color_dict[aa] = color_

		self.default_aa_size = 2
		self.master_rad = 10
		self.offset = 0
		self.mastrix = MatrixInfo.blosum80

		#Helix type-specific parameters
		self.master_angle = 102 if coiled_coil else 100
		self.step = 7 if coiled_coil else 18 #Number of aa's per wheel level
		self.width = 8 if coiled_coil else 2 #Number of levels
		self.wheel_sep = 80 if coiled_coil else 60
		self.level_sep = 2.2 if coiled_coil else 2

class Residue:
	def __init__(self, _id, _size, _num):
		self.id = _id
		self.size = _size
		self.num = _num
	
	def __str__(self):
		return self.identity

class Wheel_set:
	class Wheel:
		class Level:
			def __init__(self, _residues, _complete):
				self.residues = _residues
				self.complete = _complete

			def __str__(self):
				_ = ""
				return _.join([resi.identity for resi in self.residues])

		


		def __init__(self, rules_, identities_, scores_, wheel_num_, first_resi_num_):
			self.rules = rules_
			self.first_resi_num = first_resi_num_
			self.wheel_num = wheel_num_
			self.identities = identities_
			self.scores = scores_			
			self.size = len(self.identities)
			self.residues = self.get_residues()
			self.levels = self.get_levels()
			self.source = self.build_source()

		def __str__(self):
			return (str(self.identities) + " , wheel_num = " + str(self.wheel_num))
			

		def get_residues(self):
			_residues = []
			for i in range(self.size):
				_num = self.first_resi_num + i
				_residues.append(self.Residue(self.identities[i], self.scores[i], _num))
			# return [self.Residue(self.identities[i], self.scores[i]) for i in range(self.size)]
			return _residues
			
		def get_levels(self):
			_levels = []
			step = self.rules.step
			last_delim = 0

			for delim in range(step, self.size, step):
				_levels.append(self.Level(self.residues[last_delim:delim], True))
				last_delim = delim
			
			if last_delim != self.size:
				_levels.append(self.Level(self.residues[last_delim:], False))

			return _levels

		#Get x and y coordinates based on angle and ring number (level)
		def get_coor(self, angle, level = 0, degrees = True):
			master_rad = self.rules.master_rad
			level_sep = self.rules.level_sep

			if degrees:
				angle = math.radians(angle)
			
			rad = master_rad + (level * level_sep)

			x = rad * math.sin(angle)
			y = rad * math.cos(angle)

			return x, y

		#Build column data source to be used to construct bokeh diagram
		def build_source(self):
			colors, Xs, Ys, names, aa_num, mean_scores = ([] for _ in range(6))

			#Might want to rework system for converting conservation score to circle size
			min = abs(np.amin(self.scores))
			min = 9.4

			for num, aa in enumerate(self.identities):
				color_dict = self.rules.color_dict
				master_angle = self.rules.master_angle
				wheel_sep = self.rules.wheel_sep
				step = self.rules.step

				colors.append(color_dict[aa])
				print("step =", step)
				print("wheel_num =", self.wheel_num)
				print("num =", num)

				score_ = ((self.scores[num + (self.wheel_num * step)] + min + 2.5)  * 4)
				mean_scores.append(score_)

				angle = master_angle * num
				level = math.floor(num / step)

				x_, y_ = self.get_coor(angle, level = level)

				Xs.append(x_)

				#Distance between wheels
				Ys.append(y_ - (self.wheel_num * wheel_sep))
				names.append(aa)
				aa_num.append(num + (self.wheel_num * step) + 1)

			source = ColumnDataSource(data = dict(_Xs = Xs, _Ys = Ys, _names = names, _colors = colors, _aa_num = aa_num, _scores = mean_scores))

			input()
			return source

		def plot(self, plot):
			
		#Generate helical wheel projection
			step = self.rules.step
			data = self.source.data
			start = self.first_resi_num
			stop = start + self.size

			plot.line(x = data["_Xs"][:step], y = data["_Ys"][:step])
			plot.circle(x =  "_Xs", y =  "_Ys", color =  "_colors", name = "_names", size = "_scores", source = source)

			labels = LabelSet(
								x = "_Xs", y = "_Ys", text = "_names", level = "overlay", x_offset = -5, 
								y_offset = -5, source = self.source, render_mode = "canvas", text_color = "white"
							)

			seq_range = Label(
								x = data["_Xs"][0] + 4, y = data["_Ys"][0] - 27,
								text= ("Sequence Range: " + str(start + 1) + "-" + str(stop)),
								text_font_size = "18pt"
							)

			plot.add_layout(labels)
			plot.add_layout(seq_range)

	def __init__(self, rules_, identities_, scores_, wheels_ = None):
		self.rules = rules_
		self.identities = identities_
		self.scores = scores_

		if wheels_ is None:
			self.wheels = self.get_wheels()

	def get_wheels(self):
		_wheels = []
		bounds_set = self.get_bounds(self.identities)

		for wheel_num, bounds in enumerate(bounds_set):
			sub_ids = self.identities[bounds[0]:bounds[1]]
			sub_scores = self.scores[bounds[0]:bounds[1]]
			_wheels.append(self.Wheel(self.rules, sub_ids, sub_scores, wheel_num, bounds[0]))

		return _wheels

	#Get aa start and stop indices corresponding to each ring in HWP
	def get_bounds(self, seq):
		step = self.rules.step
		width = self.rules.width

		if (step * width) > len(seq):
			return [(0, len(seq))]

		bounds = []
		for start in range(0, len(seq), step):
			if start + (step * width) > len(seq):
				bounds.append((start, len(seq)))
			else:
				bounds.append((start, start + (step * width)))
			
		return bounds

	

	def plot_wheels(self, plot):
		for wheel in self.wheels:
			wheel.plot(plot)

#Assesses a user-defined subset of the provided sequence and alignment
def handle_range(seq, alignment):
	range_ = args.seq_range

	if range_ is None:
		return seq, alignment

	range_search = re.search("(\d+)-(\d+)", range_)

	try:
		start = int(range_search.group(1))
		stop = int(range_search.group(2))

		if start < 0 or stop > len(seq) or stop < start:
			raise Exception
	except:
		raise ValueError("Invalid range provided.")

	new_seq = seq[start:stop]

	new_align = []
	if alignment is not None:
		for sequence in alignment:
			new_align.append(sequence[start:stop])

	return new_seq, new_align

#Return conservation score from Blosum matrix for given aa pair
def get_score(pair, matrix):
	if pair in matrix:
		return matrix[pair]

	return matrix[tuple(reversed(pair))]

#Constructs matrix of conservation scores for each amino acid position in each alignment sequence
def get_score_matrix(align_seqs, matrix):
	gap_cost = -10
	ref_seq = align_seqs[0]
	score_matrix = np.zeros((len(align_seqs) - 1, len(align_seqs[0])))

	for seq_num, seq in enumerate(align_seqs[1:]):
		for aa_num in range(len(align_seqs[0])):
			pair = (ref_seq[aa_num], seq[aa_num])
			if "-" in pair:
				score_matrix[seq_num, aa_num] = gap_cost
			else:
				score_matrix[seq_num, aa_num] = get_score(pair, matrix)

	return score_matrix


#Get x and y coordinates based on angle and ring number (level)
def get_coor(angle, level = 0, degrees = True):
	global master_rad
	global level_sep

	if degrees:
		angle = math.radians(angle)
	
	rad = master_rad + (level * level_sep)

	x = rad * math.sin(angle)
	y = rad * math.cos(angle)

	return x, y

#Build column data source to be used to construct bokeh diagram
def build_CDS(seq, wheel_num, scores):
	colors, Xs, Ys, names, aa_num, mean_scores = ([] for _ in range(6))

	#Might want to rework system for converting conservation score to circle size
	min = abs(np.amin(scores))
	min = 9.4

	for num, aa in enumerate(seq):
		global color_dict
		global master_angle
		global wheel_sep

		colors.append(color_dict[aa])
		score_ = ((scores[num + (wheel_num * step)] + min + 2.5)  * 4)
		mean_scores.append(score_)

		angle = master_angle * num
		level = math.floor(num / step)

		x_, y_ = get_coor(angle, level = level)

		Xs.append(x_)

		#Distance between wheels
		Ys.append(y_ - (wheel_num * wheel_sep))
		names.append(aa)
		aa_num.append(num + (wheel_num * step) + 1)

	source = ColumnDataSource(data = dict(_Xs = Xs, _Ys = Ys, _names = names, _colors = colors, _aa_num = aa_num, _scores = mean_scores))

	return source



def main():
	blosum = MatrixInfo.blosum80
	global seq
	global align_seqs

	seq, align_seqs = handle_range(seq, align_seqs)

	hover = HoverTool(tooltips = [("aa Position", "@_aa_num"), ("conservation score", "@_scores")])

	plot = figure(
					tools = [hover, "box_select, box_zoom, wheel_zoom, pan, reset, save"],
					x_range = (-20, 20), y_range = (-20, 20),
					height = 800, width = 800
				)

	# bounds = get_bounds(seq)

	if args.align_path is None:
		#SETS DEFAULT AA SIZE
		mean_scores = np.full((len(seq)), 2)
	else:
		score_matrix = get_score_matrix(align_seqs, blosum)
		mean_scores = np.mean(score_matrix, axis = 0)

	#-----------------------------------------------------------


	master_wheel = Wheel_set(Rule_set(args.coiled_coil), seq, mean_scores)

	master_wheel.plot_wheels(plot)

	#-----------------------------------------------------------


	# for num, couple in enumerate(bounds):
	# 	source = build_CDS(seq[couple[0]:couple[1]], num, mean_scores)
	# 	plot_wheel(source, couple[0], couple[1], plot)

	plot.grid.visible = False
	plot.axis.visible = False


	show(plot)

main()

#------------TODOS--------------
"""
	--Clean up variable names
	--Move out of gloabl space
	--Make colors linear, proportional to hydrophobicity not categorical
	--Add button to snap view to next wheel
	--Add graphical user interface?
	--Handle * at the end
"""
