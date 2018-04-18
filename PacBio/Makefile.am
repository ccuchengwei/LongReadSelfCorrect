noinst_LIBRARIES = libpacbio.a

libpacbio_a_CPPFLAGS = \
	-I$(top_srcdir)/Algorithm \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/SQG \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/Thirdparty \
	-I$(top_srcdir)/Concurrency \
	-I$(top_srcdir)/FMIndexWalk \
	-I$(top_srcdir)/Thirdparty/rollinghash 

libpacbio_a_SOURCES = \
	SAIPBSelfCTree.h SAIPBSelfCTree.cpp \
	SAIPBHybridCTree.h SAIPBHybridCTree.cpp \
	PacBioSelfCorrectionProcess.h PacBioSelfCorrectionProcess.cpp \
	PacBioHybridCorrectionProcess.h PacBioHybridCorrectionProcess.cpp \
        LongReadOverlap.h LongReadOverlap.cpp \
        LongReadExtend.h LongReadExtend.cpp \
        LongReadProbe.h LongReadProbe.cpp \
	PBOverlapTree.h PBOverlapTree.cpp \
	ShortReadOverlapTree.h ShortReadOverlapTree.cpp \
	HashtableSearch.h HashtableSearch.cpp \
	LongReadCorrectByOverlap.h LongReadCorrectByOverlap.cpp \
	CheckKmerProcess.h CheckKmerProcess.cpp \
	IntervalTree.h IntervalTree.cpp IntervalTreeInstantiation.cpp \
	KmerThreshold.h KmerThreshold.cpp \
	SeedFeature.h SeedFeature.cpp \
	KmerFeature.h KmerFeature.cpp
