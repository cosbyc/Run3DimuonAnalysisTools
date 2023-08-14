#!/usr/bin/python3                                                                                                                                            
#-----------------------------------------------   

import ROOT
from ROOT import *
from math import *
#import gfal2                                                                                                      
import sys, os, pwd, subprocess, glob, fnmatch
import optparse, shlex, re
import time
from time import gmtime, strftime
import math
from array import array

Z_MASS = 91.1876
ETA_MASS = 0.547862
MU_MASS = 0.105658

#define function for parsing options
def parseOptions():

	usage = ('usage: %prog [options]\n'
			+ '%prog -h for help')
	parser = optparse.OptionParser(usage)

	# input options
	parser.add_option('-i', '--input', dest='INPUT', type='string', help='input file')
	parser.add_option('-o', '--output', dest='OUTPUT', type='string', help='output file')
	parser.add_option('-n', '--njobs', dest='NJOBS', type=int, help='njobs')
	parser.add_option('-j', '--job', dest='JOB', type=int, help='job')

	# store options and arguments as global variables
	global opt, args
	(opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quite = 0):
	status, output = subprocess.getstatusoutput(cmd)
	if (status !=0 and not quite):
		print('Error in processing command:\n 	['+cmd+']')
		print('Output:\n   ['+output+'] \n')
	return output

def fillHistogram():

	global opt, args
	parseOptions()

	ROOT.gROOT.SetBatch()
	print(opt.INPUT)

	listDir = "/afs/cern.ch/user/j/jfriesen/CMSSW_12_4_2/src/Run3DimuonAnalysisTools"	
	files = [ line for line in open(listDir+"/muMuGammaTree_ntuples.txt")]

	N = len(files)

	first = int(float(N)/float(opt.NJOBS)*float(opt.JOB-1))
	last = int(float(N)/float(opt.NJOBS)*float(opt.JOB))

	print(first, last)

	t_scoutMuon = TChain("scoutingTree/tree")
	for i in range(len(files)):
		if (i<first or i>=last): continue
		print(files[i])
		t_scoutMuon.Add(files[i])
		print(t_scoutMuon.GetEntries())

	eta_low = 0
	eta_high = 1
	eta_bins = 500

	bin_width = 0.001
	mmg_low = round(-15./bin_width)*bin_width
	mmg_high = round(15./bin_width)*bin_width
	mmg_bins = round((mmg_high - mmg_low)/bin_width)

	photon_collections = ["slimmedPhotons","pfCandPhotons","pfCandPhotonsPtMax10","pfCandPhotonsPtMin10"]
	selections = ["closestToEta","minDr","minDrEt2"]
	plots = ["massMMG","massDimu","massDiff"]
	cuts = ["all","isEta","isNotEta","massDimuMax0p528"]

	config = {}
	for collection in photon_collections :
		config[collection] = {}
		for selection in selections :
			config[collection][selection] = {}
			for plot in plots :
				config[collection][selection][plot] = {}
				for cut in cuts :
					config[collection][selection][plot][cut] = ROOT.TH1F(collection+"_"+selection+"_"+plot+"_"+cut,collection+"_"+selection+"_"+plot+"_"+cut,mmg_bins,mmg_low,mmg_high)
			config[collection][selection]["best value"] = -1
			config[collection][selection]["best photon"] = -1

	rand = ROOT.TRandom()

	verbose = False
	anyMassMMG = True

	i_event = 0
	for ev in t_scoutMuon :

		if(i%10000==0): print("event:",i_event)
		if (verbose) : print("\nevent:",i_event)
		i_event+=1

		mu1 = ROOT.Math.PtEtaPhiMVector(ev.pt1, ev.eta1, ev.phi1, MU_MASS) 
		mu2 = ROOT.Math.PtEtaPhiMVector(ev.pt2, ev.eta2, ev.phi2, MU_MASS)
		if (verbose) : print("mu1.M",mu1.M(),"mu2.M",mu2.M())
		dimu = mu1+mu2
		mass_dimu = dimu.M()
		if (verbose) : print("mass_dimu",mass_dimu)

		if (verbose) : print("slimmedPhotons")
		for i_photon in range(len(ev.slimmedPhotonPt)) :
			#print("	photon",i_photon)
			collection = "slimmedPhotons"
			gamma = ROOT.Math.PtEtaPhiMVector(ev.slimmedPhotonPt[i_photon], ev.slimmedPhotonEta[i_photon], ev.slimmedPhotonPhi[i_photon], 0)
			mass_mmg = (dimu + gamma).M()
			dr = abs(ROOT.Math.VectorUtil.DeltaR(dimu, gamma))
			if (verbose) : print("	photon massMMG", mass_mmg)
			for selection in selections :
				if (verbose) : print("		selection",selection)
				if selection == "closestToEta" : value = abs(mass_mmg - ETA_MASS)
				if selection == "minDr" : value = dr
				if selection == "minDrEt2" : value = dr / (ev.slimmedPhotonPt[i_photon] ** 2)
				if config[collection][selection]["best value"] == -1 or value < config[collection][selection]["best value"] :
					if (verbose) : print("			best value",value,i_photon,collection,selection)
					config[collection][selection]["best value"] = value
					config[collection][selection]["best photon"] = i_photon

		if (verbose) : print("pfCandPhotons")
		for i_photon in range(len(ev.photonPt)) :
			#print("	photon",i_photon)
			gamma = ROOT.Math.PtEtaPhiMVector(ev.photonPt[i_photon], ev.photonEta[i_photon], ev.photonPhi[i_photon], 0)
			mass_mmg = (dimu + gamma).M()
			dr = abs(ROOT.Math.VectorUtil.DeltaR(dimu, gamma))
			if (verbose) : print("	photon massMMG", mass_mmg)
			for collection in photon_collections[1:] :
				if (verbose) : print("		collection", collection)
				if collection == "pfCandPhotonsPtMax10" and ev.photonPt[i_photon] > 10 : continue
				if collection == "pfCandPhotonsPtMin10" and ev.photonPt[i_photon] < 10 : continue
				for selection in selections :
					if (verbose) : print("			selection", selection)
					if selection == "closestToEta" : value = abs(mass_mmg - ETA_MASS)
					if selection == "minDr" : value = dr
					if selection == "minDrEt2" : value = dr / ev.photonEt2[i_photon]
					if config[collection][selection]["best value"] == -1 or value < config[collection][selection]["best value"] :
						if (verbose) : print("				best value",value,i_photon,collection,selection)
						config[collection][selection]["best value"] = value
						config[collection][selection]["best photon"] = i_photon


		for collection in photon_collections : 
			for selection in selections :
				if verbose : print(collection,selection,config[collection][selection]["best value"],config[collection][selection]["best photon"])
				if config[collection][selection]["best photon"] < 0 : continue
				if collection == "slimmedPhotons" :
					gamma = ROOT.Math.PtEtaPhiMVector(ev.slimmedPhotonPt[config[collection][selection]["best photon"]], ev.slimmedPhotonEta[config[collection][selection]["best photon"]], ev.slimmedPhotonPhi[config[collection][selection]["best photon"]], 0)
				else : 
					gamma = ROOT.Math.PtEtaPhiMVector(ev.photonPt[config[collection][selection]["best photon"]], ev.photonEta[config[collection][selection]["best photon"]], ev.photonPhi[config[collection][selection]["best photon"]], 0)
				config[collection][selection]["best photon"] = -1
				config[collection][selection]["best value"] = -1
				mass_mmg = (dimu+gamma).M()
				for cut in cuts :
					if cut == "isEta" and abs(mass_mmg - ETA_MASS) > abs(mass_dimu - ETA_MASS) : continue
					if cut == "isNotEta" and abs(mass_mmg - ETA_MASS) < abs(mass_dimu - ETA_MASS) : continue
					if cut == "massDimuMax0p528" and mass_dimu > 0.528 : continue
					config[collection][selection]["massMMG"][cut].Fill(mass_mmg)
					if(mass_mmg > 0.5 and mass_mmg < 0.6) :
						config[collection][selection]["massDimu"][cut].Fill(mass_dimu)
						config[collection][selection]["massDiff"][cut].Fill(abs(mass_mmg - ETA_MASS) - abs(mass_dimu - ETA_MASS))


	print("saving as "+str(opt.OUTPUT)+str(opt.JOB)+".root")
	outfile = ROOT.TFile(str(opt.OUTPUT)+str(opt.JOB)+".root", "recreate")

	for collection in photon_collections :
		for selection in selections :
			for plot in plots :
				for cut in cuts :
					print("saving as",collection+"_"+selection+"_"+plot+"_"+cut,"with",(config[collection][selection][plot][cut]).GetEntries(),"entries")
					outfile.WriteObject(config[collection][selection][plot][cut], collection+"_"+selection+"_"+plot+"_"+cut)


	outfile.Close()

if __name__ == "__main__":
	fillHistogram()


