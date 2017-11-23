///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess.cpp - Self-correction using FM-index walk for PacBio reads
//

#ifndef PACBIOSELFCORRECTION_H
#define PACBIOSELFCORRECTION_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "BWTAlgorithms.h"

// functions

//
int PacBioSelfCorrectionMain(int argc, char** argv);

// options
void parsePacBioSelfCorrectionOptions(int argc, char** argv);

#endif
