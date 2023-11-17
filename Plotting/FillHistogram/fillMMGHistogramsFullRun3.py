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
	parser.add_option('-l', dest='LIST', help='ntuple list', default="/afs/cern.ch/user/j/jfriesen/CMSSW_13_0_10/src/Run3DimuonAnalysisTools/Plotting/FillHistogram/muMuGammaTree_ntuples_fullRun3.txt")

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

	files = [ "root://cmsxrootd.fnal.gov//"+line for line in open(str(opt.LIST))]

	N = len(files)

	first = int(float(N)/float(opt.NJOBS)*float(opt.JOB-1))
	last = int(float(N)/float(opt.NJOBS)*float(opt.JOB))

	print(first, last)

	t_scoutMuon = TChain("tree/tree")
	for i in range(len(files)):
		if (i<first or i>=last): continue
		print(files[i])
		t_scoutMuon.Add(files[i])
		print(t_scoutMuon.GetEntries())

	eta_low = 0
	eta_high = 1
	eta_bins = 500

	bin_width = 0.001
	mmg_low = round(0./bin_width)*bin_width
	mmg_high = round(10./bin_width)*bin_width
	mmg_bins = round((mmg_high - mmg_low)/bin_width)

	photon_collections = ["slimmedPhotons","pfCandPhotons","slimmedOrPfCandPhotons"]
	selections = ["minDr"]
	plots = ["massMMG","massDimu","massDimu_massMMG0p503to0p571"]
	cuts = ["all","isEta"]

	verbose = False
	muon_selection = True

	dimuMass = ROOT.TH1F("massDimu","massDimu",mmg_bins,mmg_low,mmg_high)

	config = {}
	for collection in photon_collections :
		config[collection] = {}
		for selection in selections :
			config[collection][selection] = {}
			for plot in plots :
				config[collection][selection][plot] = {}
				for cut in cuts :
					if verbose : print("	"+collection+"_"+selection+"_"+plot+"_"+cut)
					config[collection][selection][plot][cut] = ROOT.TH1F(collection+"_"+selection+"_"+plot+"_"+cut,collection+"_"+selection+"_"+plot+"_"+cut,mmg_bins,mmg_low,mmg_high)
			config[collection][selection]["best value"] = -1
			config[collection][selection]["best photon"] = -1

	rand = ROOT.TRandom()

	i_event = 0
	for ev in t_scoutMuon :

		if(i_event%10000==0): print("event:",i_event)
		if (verbose) : print("\nevent:",i_event)
		i_event+=1

		if (muon_selection) and (not ev.muonID1[3] or not ev.muonID2[3] or ev.probVtx < 0.05) : continue

		mu1 = ROOT.Math.PtEtaPhiMVector(ev.pt1, ev.eta1, ev.phi1, MU_MASS) 
		mu2 = ROOT.Math.PtEtaPhiMVector(ev.pt2, ev.eta2, ev.phi2, MU_MASS)
		if (verbose) : print("mu1.M",mu1.M(),"mu2.M",mu2.M())
		dimu = mu1+mu2
		mass_dimu = dimu.M()
		if (verbose) : print("mass_dimu",mass_dimu)
		dimuMass.Fill(mass_dimu)

		slimmedOrPfCandPhotons = ""

		if ("slimmedOrPfCandPhotons" in photon_collections) or ("pfCandPhotons" in photon_collections) :
			if (verbose) : print("pfCandPhotons")
			for i_photon in range(len(ev.pfCandPhotonPt)) :
				#print("	photon",i_photon)
				gamma = ROOT.Math.PtEtaPhiMVector(ev.pfCandPhotonPt[i_photon], ev.pfCandPhotonEta[i_photon], ev.pfCandPhotonPhi[i_photon], 0)
				mass_mmg = (dimu + gamma).M()
				dr = abs(ROOT.Math.VectorUtil.DeltaR(dimu, gamma))
				if (verbose) : print("	photon massMMG", mass_mmg)
				for collection in photon_collections :
					if (collection == "slimmedPhotons"): continue
					if (verbose) : print("		collection", collection)
					if collection == "pfCandPhotonsPtMax10" and ev.pfCandPhotonPt[i_photon] > 10 : continue
					if collection == "pfCandPhotonsPtMin10" and ev.pfCandPhotonPt[i_photon] < 10 : continue
					for selection in selections :
						if selection == "closestToEta" : value = abs(mass_mmg - ETA_MASS)
						if selection == "minDr" : value = dr
						if selection == "minDrEt2" : value = dr / ev.pfCandPhotonEt2[i_photon]
						if (verbose) : print("			selection", selection, "value", value, "best value", config[collection][selection]["best value"])
						if config[collection][selection]["best value"] == -1 or value < config[collection][selection]["best value"] :
							if (verbose) : print("				best value",value,i_photon,collection,selection)
							config[collection][selection]["best value"] = value
							config[collection][selection]["best photon"] = i_photon
							if collection == "slimmedOrPfCandPhotons" : slimmedOrPfCandPhotons = "pfCandPhotons"
		
		if ("slimmedPhotons" in photon_collections) or ("slimmedOrPfCandPhotons" in photon_collections) :
			if (verbose) : print("slimmedPhotons")
			for i_photon in range(len(ev.slimmedPhotonPt)) :
				#print("	photon",i_photon)
				collection = "slimmedPhotons"
				gamma = ROOT.Math.PtEtaPhiMVector(ev.slimmedPhotonPt[i_photon], ev.slimmedPhotonEta[i_photon], ev.slimmedPhotonPhi[i_photon], 0)
				mass_mmg = (dimu + gamma).M()
				dr = abs(ROOT.Math.VectorUtil.DeltaR(dimu, gamma))
				if (verbose) : print("	photon massMMG", mass_mmg)
				for selection in selections :
					if (verbose) : print("		selection", selection, "value", value)
					if selection == "closestToEta" : value = abs(mass_mmg - ETA_MASS)
					if selection == "minDr" : value = dr
					if selection == "minDrEt2" : value = dr / (ev.slimmedPhotonPt[i_photon] ** 2)
					if config[collection][selection]["best value"] == -1 or value < config[collection][selection]["best value"] :
						if (verbose) : print("			best value",value,i_photon,collection,selection)
						config[collection][selection]["best value"] = value
						config[collection][selection]["best photon"] = i_photon
			if len(ev.slimmedPhotonPt) > 0 and ("slimmedOrPfCandPhotons" in photon_collections) :
				if (verbose) : print("slimmedOrPfCandPhotons from slimmed")
				for selection in selections :
					if (config["slimmedOrPfCandPhotons"][selection]["best value"] == -1) or (config["slimmedPhotons"][selection]["best value"] < config["slimmedOrPfCandPhotons"][selection]["best value"]) :
						if (verbose) : print("				best value",value,i_photon,collection,selection)
						config["slimmedOrPfCandPhotons"][selection]["best value"] = config["slimmedPhotons"][selection]["best value"]
						config["slimmedOrPfCandPhotons"][selection]["best photon"] = config["slimmedPhotons"][selection]["best photon"]
						slimmedOrPfCandPhotons = "slimmedPhotons"

		for collection in photon_collections : 
			for selection in selections :
				if verbose : print(collection,selection,config[collection][selection]["best value"],config[collection][selection]["best photon"])
				if config[collection][selection]["best photon"] < 0 : continue
				if collection == "slimmedPhotons" or (collection == "slimmedOrPfCandPhotons" and slimmedOrPfCandPhotons == "slimmedPhotons"):
					gamma = ROOT.Math.PtEtaPhiMVector(ev.slimmedPhotonPt[config[collection][selection]["best photon"]], ev.slimmedPhotonEta[config[collection][selection]["best photon"]], ev.slimmedPhotonPhi[config[collection][selection]["best photon"]], 0)
				else : 
					gamma = ROOT.Math.PtEtaPhiMVector(ev.pfCandPhotonPt[config[collection][selection]["best photon"]], ev.pfCandPhotonEta[config[collection][selection]["best photon"]], ev.pfCandPhotonPhi[config[collection][selection]["best photon"]], 0)
				config[collection][selection]["best photon"] = -1
				config[collection][selection]["best value"] = -1
				mass_mmg = (dimu+gamma).M()
				for cut in cuts :
					if cut == "isEta" and abs(mass_mmg - ETA_MASS) > abs(mass_dimu - ETA_MASS) : continue
					config[collection][selection]["massMMG"][cut].Fill(mass_mmg)
					config[collection][selection]["massDimu"][cut].Fill(mass_dimu)
					if(mass_mmg > 0.503 and mass_mmg < 0.571) :
						config[collection][selection]["massDimu_massMMG0p503to0p571"][cut].Fill(mass_dimu)


	print("saving as "+str(opt.OUTPUT)+str(opt.JOB)+".root")
	outfile = ROOT.TFile(str(opt.OUTPUT)+str(opt.JOB)+".root", "recreate")

	print("saving as massDimu with",dimuMass.GetEntries(),"entries")
	outfile.WriteObject(dimuMass, "massDimu")

	for collection in photon_collections :
		for selection in selections :
			for plot in plots :
				for cut in cuts :
					print("saving as",collection+"_"+selection+"_"+plot+"_"+cut,"with",(config[collection][selection][plot][cut]).GetEntries(),"entries")
					outfile.WriteObject(config[collection][selection][plot][cut], collection+"_"+selection+"_"+plot+"_"+cut)


	outfile.Close()

if __name__ == "__main__":
	fillHistogram()