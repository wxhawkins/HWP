import math
from functools import partial
import numpy as np
import re

from bokeh.models import HoverTool, Button, Label, LabelSet
from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh import palettes

from Bio.SubsMat import MatrixInfo
from Bio import SeqIO

import argparse
import time

#Establish argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-s", action = "store", type = str, dest = "seq_path", required = True)
parser.add_argument("-a", action = "store", type = str, dest = "align_path", default = None)
parser.add_argument("-r", action = "store", type = str, dest = "seq_range", default = None)
parser.add_argument("-c", action = "store_true", dest = "coiled_coil", default = False)
args = parser.parse_args()

#Master variables
if args.coiled_coil:
	master_rad = 10
	master_angle = 102
	step = 7 #Number of aa's per wheel level
	width = 8 #Number of levels
	blosum = MatrixInfo.blosum80
	wheel_sep = 80
	level_sep = 2.2
else:
	master_rad = 10
	master_angle = 100
	step = 18 #Number of aa's per wheel level
	width = 2 #Number of levels
	blosum = MatrixInfo.blosum80
	wheel_sep = 60
	level_sep = 2

#Establish alphabet
all_aas = "WIFLCMVYPATHGSQNEDKR"

#Set colors for each amino acid based on hydrophobicity (two aa's per color)
pal = palettes.all_palettes['RdYlBu'][10]
color_dict = {}
for i, aa in enumerate(all_aas):
	color_ = pal[math.floor(i / 2)]
	color_dict[aa] = color_

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

#Assesses a user-defined subset of the provieded sequence and alignment
def handle_range(seq, alignment):
	range = args.seq_range

	if range is None:
		return seq, alignment

	range_search = re.search("(\d+)-(\d+)", range)

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
		print("step =", step)
		print("wheel_num =", wheel_num)
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

	input()
	return source

#Generate helical wheel projection
def plot_wheel(source, start, stop, plot):
	global step
	data = source.data

	plot.line(x = data["_Xs"][:step], y = data["_Ys"][:step])
	plot.circle(x =  "_Xs", y =  "_Ys", color =  "_colors", name = "_names", size = "_scores", source = source)

	labels = LabelSet(
						x = "_Xs", y = "_Ys", text = "_names", level = "overlay", x_offset = -5, 
						y_offset = -5, source = source, render_mode = "canvas", text_color = "white"
					 )

	seq_range = Label(
						x = data["_Xs"][0] + 4, y = data["_Ys"][0] - 27,
						text= ("Sequence Range: " + str(start + 1) + "-" + str(stop)),
						text_font_size = "18pt"
					 )

	plot.add_layout(labels)
	plot.add_layout(seq_range)

#Get aa start and stop indices corresponding to each ring in HWP
def get_bounds(seq):
	global step
	global width

	if (step * width) > len(seq):
		return [(0, len(seq))]

	bounds = []

	for start in range(0, len(seq), step):
		if start + (step * width) > len(seq):
			bounds.append((start, len(seq)))
		else:
			bounds.append((start, start + (step * width)))
		
	return bounds

def main():
	global seq
	global align_seqs

	seq, align_seqs = handle_range(seq, align_seqs)

	hover = HoverTool(tooltips = [("aa Position", "@_aa_num"), ("conservation score", "@_scores")])

	plot = figure(
					tools = [hover, "box_select, box_zoom, wheel_zoom, pan, reset, save"],
					x_range = (-20, 20), y_range = (-20, 20),
					height = 800, width = 800
				)

	bounds = get_bounds(seq)

	if args.align_path is None:
		#SETS DEFAULT AA SIZE
		mean_scores = np.full((len(seq)), 2)
	else:
		score_matrix = get_score_matrix(align_seqs, blosum)
		mean_scores = np.mean(score_matrix, axis = 0)

	# for i, score in enumerate(mean_scores):
	# 	print(i+1, " = ", score)

	for num, couple in enumerate(bounds):
		source = build_CDS(seq[couple[0]:couple[1]], num, mean_scores)
		plot_wheel(source, couple[0], couple[1], plot)

	plot.grid.visible = False
	plot.axis.visible = False

	# print(palettes.all_palettes['RdYlBu'][10])

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
