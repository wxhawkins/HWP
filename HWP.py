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

# Establish argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-s", action = "store", type = str, dest = "seq_path", required = True)
parser.add_argument("-a", action = "store", type = str, dest = "align_path", default = None)
parser.add_argument("-r", action = "store", type = str, dest = "seq_range", default = None)
parser.add_argument("-c", action = "store_true", dest = "coiled_coil", default = False)
parser.add_argument("-f", action = "store_true", dest = "find", default = False)
args = parser.parse_args()

class Rule_set:
	def __init__(self, coiled_coil = False):
		# Establish amino acid alphabet
		self.all_aas = "WIFLCMVYPATHGSQNEDKR"

		self.hydrophobe_dict = {
								"A":0.310, "D":-0.770, "E":-0.640, "I":1.800, "M":1.230, "S":-0.040, "Y": 0.960,
								"R":-1.010, "C":1.540, "G":0.000, "L":1.700, "F":1.790, "T":0.260, "V":1.200,
								"N":-0.600, "Q":-0.220, "H":0.130, "K":-0.990, "P":0.720, "W":2.250
						  	   }

		# Order of amino acids when viewed as helical wheel
		self.wheel_order = [0, 11, 4, 15, 8, 1, 12, 5, 16, 9, 2, 13, 6, 17, 10, 3, 14, 7]

		# Set colors for each amino acid based on hydrophobicity (two aa's per color)
		pal = palettes.all_palettes['RdYlBu'][10]
		self.color_dict = {}
		for i, aa in enumerate(self.all_aas):
			color_ = pal[math.floor(i / 2)]
			self.color_dict[aa] = color_

		self.default_resi_size = 2
		self.master_rad = 10
		self.offset = 0
		self.matrix = MatrixInfo.blosum80

		# Helix type-specific parameters
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
		return self.id

class Wheel_set:
	def __init__(self, _rules, _residues, _wheels=None):
		self.rules = _rules
		self.residues = _residues
		self.length = len(self.residues)
		self.wheels = _wheels

		if _wheels is None:
			self.wheels = self.get_wheels()

	def get_wheels(self):
		_wheels = []
		bounds_set = self.get_bounds()

		for wheel_num, bounds in enumerate(bounds_set):
			_wheels.append(Wheel(self.rules, self.residues[bounds[0]:bounds[1]], wheel_num))

		return _wheels

	# Get aa start and stop indices corresponding to each ring in HWP
	def get_bounds(self):
		step = self.rules.step
		width = self.rules.width
		length = self.length

		if (step * width) > length:
			return [(0, length)]

		bounds = []
		for start in range(0, length, step):
			if start + (step * width) > length:
				bounds.append((start, length))
			else:
				bounds.append((start, start + (step * width)))
			
		return bounds

	def plot_wheels(self, plot):
		for wheel in self.wheels:
			wheel.plot(plot)

class Wheel:
	class Level:
		def __init__(self, _residues, _complete):
			self.residues = _residues
			self.complete = _complete

		def __str__(self):
			_ = ""
			return _.join([resi.id for resi in self.residues])

	def __init__(self, _rules, _residues, _wheel_num):
		self.rules = _rules
		self.wheel_num = _wheel_num
		self.residues = _residues
		self.length = len(self.residues)
		self.levels = self.get_levels()
		self.source = self.build_source()

	def __str__(self):
		aa_string = ""
		aa_string = aa_string.join([resi.id for resi in self.residues])
		return (aa_string + " , wheel_num = " + str(self.wheel_num))
		
	def get_levels(self):
		_levels = []
		step = self.rules.step
		last_delim = 0

		for delim in range(step, self.length, step):
			_levels.append(self.Level(self.residues[last_delim:delim], True))
			last_delim = delim
		
		if last_delim != self.length:
			_levels.append(self.Level(self.residues[last_delim:], False))

		return _levels

	# Get x and y coordinates based on angle and ring number (level)
	def get_coor(self, angle, level = 0, degrees = True):
		master_rad = self.rules.master_rad
		level_sep = self.rules.level_sep

		if degrees:
			angle = math.radians(angle)
		
		rad = master_rad + (level * level_sep)

		x = rad * math.sin(angle)
		y = rad * math.cos(angle)

		return x, y

	# Build column data source to be used to construct bokeh diagram
	def build_source(self):
		colors, Xs, Ys, names, aa_num, mean_scores = ([] for _ in range(6))
		color_dict = self.rules.color_dict
		master_angle = self.rules.master_angle
		wheel_sep = self.rules.wheel_sep
		step = self.rules.step

		# Might want to rework system for converting conservation score to circle size
		scores = [resi.size for resi in self.residues]
		min = abs(np.amin(scores))
		min = 9.4

		for pos, resi in enumerate(self.residues):
			colors.append(color_dict[resi.id])

			_score = ((resi.size + min + 2.5) * 4)
			mean_scores.append(_score)

			angle = master_angle * pos
			level = math.floor(pos / step)

			x_, y_ = self.get_coor(angle, level = level)

			Xs.append(x_)

			# Distance between wheels
			Ys.append(y_ - (self.wheel_num * wheel_sep))
			names.append(resi.id)
			aa_num.append(resi.num)

		source = ColumnDataSource(data = dict(_Xs = Xs, _Ys = Ys, _names = names, _colors = colors, _aa_num = aa_num, _scores = mean_scores))

		return source

	# Generate helical wheel projection
	def plot(self, plot):
		step = self.rules.step
		data = self.source.data
		start = self.residues[0].num
		stop = self.residues[-1].num

		plot.line(x = data["_Xs"][:step], y = data["_Ys"][:step])
		plot.circle(x =  "_Xs", y =  "_Ys", color =  "_colors", name = "_names", size = "_scores", source = self.source)

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

# Gets amino acid sequence as viewed in helical wheel
def get_wheel_seq(residues):
	wheel_seq = []
	threshold = 18
	i = 0

	while i < (len(residues)):
		if len(wheel_seq) > threshold:
			wheel_seq = wheel_seq[:threshold]
			i = threshold
			threshold += 18

			if i >= len(residues):
				break

		wheel_seq.append(residues[i])
		try:
			wheel_seq.append(residues[i + 11])
		except: pass

		try:
			wheel_seq.append(residues[i + 4])
		except: pass

		try:
			wheel_seq.append(residues[i + 15])
		except: pass

		try:
			wheel_seq.append(residues[i + 8])
		except: pass
		
		i += 1

	return wheel_seq


def plot_amph_wheels(rules, residues, hp_scores, plot, wheel_seq):
	thresh = 7
	wheels = []
	for wheel_pos, score in enumerate(hp_scores):
		if score > thresh:
			if wheel_pos > 15 and wheel_pos < (len(residues) - 15): #REFINE
				wheels.append(Wheel(rules, residues[wheel_seq[wheel_pos].num-9:wheel_seq[wheel_pos].num+9], len(wheels)))
			
	for wheel in wheels:
		wheel.plot(plot)

# Assesses a user-defined subset of the provided sequence and alignment
def handle_range(residues):
	range_ = args.seq_range

	if range_ is None:
		return residues

	range_search = re.search("(\d+)-(\d+)", range_)
	
	try:
		start = int(range_search.group(1))
		stop = int(range_search.group(2))

		if start < 0 or stop > len(residues) or stop < start:
			raise Exception
	except:
		raise ValueError("Invalid range provided.")

	return residues[start:stop + 1]

# Return conservation score from Blosum matrix for given aa pair
def get_score(pair, matrix):
	if pair in matrix:
		return matrix[pair]

	return matrix[tuple(reversed(pair))]

# Constructs matrix of conservation scores for each amino acid position in each alignment sequence
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

def get_resis(seq_path, alignment_path, rules):
	# Get primary sequence to be analyzed
	try:
		seq = SeqIO.read(seq_path, "fasta").seq
		if len(seq) == 0:
			raise Exception
	except:
		raise ValueError("Invalid protein sequence provided for analyzation.")

	# Get and clean sequence alignment
	if alignment_path is None:
		mean_scores = np.full(len(seq), rules.default_resi_size)
	else:
		try:
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

			score_matrix = get_score_matrix(align_seqs, rules.matrix)
			mean_scores = np.mean(score_matrix, axis = 0)
		except:
			raise ValueError("Invalid alignment sequence provided.")

	_residues = []
	for i in range(len(seq)):
		_num = i + 1
		_residues.append(Residue(seq[i], mean_scores[i], _num))
	
	return handle_range(_residues)

def main():
	hover = HoverTool(tooltips = [("aa Position", "@_aa_num"), ("conservation score", "@_scores")])

	plot = figure(
					tools = [hover, "box_select, box_zoom, wheel_zoom, pan, reset, save"],
					x_range = (-20, 20), y_range = (-20, 20),
					height = 800, width = 800
				 )


	#-----------------------------------------------------------

	rules = Rule_set(args.coiled_coil)
	master_res_list = get_resis(args.seq_path, args.align_path, rules)

	master_wheel = Wheel_set(rules, master_res_list)
	master_wheel.plot_wheels(plot)

	#-----------------------------------------------------------
	#TEST

	# prim_seq = [resi.id for resi in master_res_list]
	# # wheel_seq = get_wheel_seq(master_res_list)

	# amph_wheels = find_phobe_patch(rules, master_res_list)

	# for wheel in amph_wheels:
	# 	wheel.plot(plot)

	# print(scores)
	# plot_amph_wheels(rules, master_res_list, scores, plot, wheel_seq)

	#-----------------------------------------------------------

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
