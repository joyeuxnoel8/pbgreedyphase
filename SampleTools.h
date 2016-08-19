#ifndef SAMPLE_TOOLS_H_
#define SAMPLE_TOOLS_H_

#include "seqan/vcf_io.h"
#include "seqan/seq_io.h"
#include "seqan/bam_io.h"
#include "seqan/align.h"
#include "seqan/basic.h"
using namespace seqan;
int LookupPhasedSampleIndex(VcfIOContext<>::TNameStore &sampleNames, string phasedSample, int &sampleIndex) {
	int i;
	for (i = 0; i < length(sampleNames); i++) {
		if (sampleNames[i] == phasedSample) {
			sampleIndex = i;
			return true;
		}
	}
	return false;
}


#endif
