//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// bwt2fa - Transform a bwt back into a set of sequences
//
#ifndef CHECKKMER_H
#define CHECKKMER_H
#include <getopt.h>
#include "config.h"

int checkkmerMain(int argc, char** argv);
void parseCHECKKMEROptions(int argc, char** argv);

#endif
