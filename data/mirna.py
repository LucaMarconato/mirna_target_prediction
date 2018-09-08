#!/usr/bin/env python3.7

class Mirna:
	def __init__(self, mirna_family):
		self.mirna_family = mirna_family

	def __hash__(self):
		return hash(self.mirna_family)

	def __eq__(self, other):
		return (self.mirna_family == other.mirna_family)

	def __ne__(self, other):
		return not(self == other)
