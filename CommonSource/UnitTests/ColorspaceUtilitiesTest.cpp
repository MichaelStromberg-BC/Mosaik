// ***************************************************************************
// ColorspaceUtilitiesTest.cpp - provides unit tests for CColorspaceUtilities.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <string>
#include "Alignment.h"
#include "ColorspaceUtilities.h"
#include "WinUnit.h"

using namespace std;

BEGIN_TEST(CColorspaceUtilities_ConvertAlignmentToBasespace) {

	// create a new instance of our colorspace utilities
	CColorspaceUtilities csu;

	// assign the reference sequence to our instance
	string bsRefSeq = "AATCTGAGAGGCAGAGGTTGCAGTGAGCTGAGATTGTGCC";
	const unsigned int bsRefSeqLen = bsRefSeq.size();

	char** pBsRefSeqs = new char*[1];
	pBsRefSeqs[0] = new char[bsRefSeqLen + 1];
	memcpy(pBsRefSeqs[0], bsRefSeq.data(), bsRefSeqLen);
	pBsRefSeqs[0][bsRefSeqLen] = 0;

	csu.SetReferenceSequences(pBsRefSeqs);

	// create some test alignments
	Alignment exact, oneSNP, oneError, snpAndError, errorAndSnp, longInsertion, shortInsertion, longDeletion, shortDeletion, shortDeletionAndError;

	exact.Reference       = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	exact.Query           = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	exact.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	exact.ReferenceBegin  = 5;
	exact.ReferenceEnd    = 33;
	exact.QueryBegin      = 0;
	exact.QueryEnd        = 28;
	exact.IsReverseStrand = false;
	exact.ReferenceIndex  = 0;

	oneSNP.Reference       = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	oneSNP.Query           = "GGGGATCGGGACACTCATCGGTGCGGGTA";
	oneSNP.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	oneSNP.ReferenceBegin  = 5;
	oneSNP.ReferenceEnd    = 33;
	oneSNP.QueryBegin      = 0;
	oneSNP.QueryEnd        = 28;
	oneSNP.IsReverseStrand = false;
	oneSNP.ReferenceIndex  = 0;

	oneError.Reference       = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	oneError.Query           = "GGGGATCGGGACACTCACCGGTGCGGGTA";
	oneError.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	oneError.ReferenceBegin  = 5;
	oneError.ReferenceEnd    = 33;
	oneError.QueryBegin      = 0;
	oneError.QueryEnd        = 28;
	oneError.IsReverseStrand = false;
	oneError.ReferenceIndex  = 0;

	snpAndError.Reference       = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	snpAndError.Query           = "GGGGATCGGGACACTCATGGGTGCGGGTA";
	snpAndError.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	snpAndError.ReferenceBegin  = 5;
	snpAndError.ReferenceEnd    = 33;
	snpAndError.QueryBegin      = 0;
	snpAndError.QueryEnd        = 28;
	snpAndError.IsReverseStrand = false;
	snpAndError.ReferenceIndex  = 0;

	errorAndSnp.Reference       = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	errorAndSnp.Query           = "GGGGATCGGGACACTGATCGGTGCGGGTA";
	errorAndSnp.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	errorAndSnp.ReferenceBegin  = 5;
	errorAndSnp.ReferenceEnd    = 33;
	errorAndSnp.QueryBegin      = 0;
	errorAndSnp.QueryEnd        = 28;
	errorAndSnp.IsReverseStrand = false;
	errorAndSnp.ReferenceIndex  = 0;

	shortInsertion.Reference       = "GGGGATCGGGACACTCG-CCGGTGCGGGTA";
	shortInsertion.Query           = "GGGGATCGGGACACTCGTGCGGTGCGGGTA";
	shortInsertion.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcd";
	shortInsertion.ReferenceBegin  = 5;
	shortInsertion.ReferenceEnd    = 33;
	shortInsertion.QueryBegin      = 0;
	shortInsertion.QueryEnd        = 29;
	shortInsertion.IsReverseStrand = false;
	shortInsertion.ReferenceIndex  = 0;

	longInsertion.Reference       = "GGGGATCGGGACACTCG----CCGGTGCGGGTA";
	longInsertion.Query           = "GGGGATCGGGACACTCGCTCTCCGGTGCGGGTA";
	longInsertion.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg";
	longInsertion.ReferenceBegin  = 5;
	longInsertion.ReferenceEnd    = 33;
	longInsertion.QueryBegin      = 0;
	longInsertion.QueryEnd        = 32;
	longInsertion.IsReverseStrand = false;
	longInsertion.ReferenceIndex  = 0;

	longDeletion.Reference       = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	longDeletion.Query           = "GGGGATCGGGACACTC----GTGCGGGTA";
	longDeletion.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXY";
	longDeletion.ReferenceBegin  = 5;
	longDeletion.ReferenceEnd    = 33;
	longDeletion.QueryBegin      = 0;
	longDeletion.QueryEnd        = 24;
	longDeletion.IsReverseStrand = false;
	longDeletion.ReferenceIndex  = 0;

	shortDeletion.Reference       = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	shortDeletion.Query           = "GGGGATCGGGACACTC-TCGGTGCGGGTA";
	shortDeletion.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZab";
	shortDeletion.ReferenceBegin  = 5;
	shortDeletion.ReferenceEnd    = 33;
	shortDeletion.QueryBegin      = 0;
	shortDeletion.QueryEnd        = 24;
	shortDeletion.IsReverseStrand = false;
	shortDeletion.ReferenceIndex  = 0;

	shortDeletionAndError.Reference       = "GGGGATCGGGACACTCGCCGGTGCGGGTA";
	shortDeletionAndError.Query           = "GGGGATCGGGACACTC-TGGGTGCGGGTA";
	shortDeletionAndError.BaseQualities   = "ABCDEFGHIJKLMNOPQRSTUVWXYZab";
	shortDeletionAndError.ReferenceBegin  = 5;
	shortDeletionAndError.ReferenceEnd    = 33;
	shortDeletionAndError.QueryBegin      = 0;
	shortDeletionAndError.QueryEnd        = 24;
	shortDeletionAndError.IsReverseStrand = false;
	shortDeletionAndError.ReferenceIndex  = 0;

	// convert the alignments
	csu.ConvertAlignmentToBasespace(exact);
	csu.ConvertAlignmentToBasespace(oneSNP);
	csu.ConvertAlignmentToBasespace(oneError);
	csu.ConvertAlignmentToBasespace(snpAndError);
	csu.ConvertAlignmentToBasespace(errorAndSnp);
	csu.ConvertAlignmentToBasespace(shortInsertion);
	csu.ConvertAlignmentToBasespace(longInsertion);
	csu.ConvertAlignmentToBasespace(longDeletion);
	csu.ConvertAlignmentToBasespace(shortDeletion);
	csu.ConvertAlignmentToBasespace(shortDeletionAndError);

	// check the results
	CMosaikString expectedReference = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	CMosaikString expectedQuery     = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	CMosaikString expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	unsigned short expectedNumMismatches = 0;
	WIN_ASSERT_EQUAL(exact.Reference, expectedReference, _T("Failed basespace conversion of the exact reference.\n"));
	WIN_ASSERT_EQUAL(exact.Query, expectedQuery, _T("Failed basespace conversion of the exact query.\n"));
	WIN_ASSERT_EQUAL(exact.BaseQualities, expectedQualities, _T("Failed basespace conversion of the exact base qualities.\n"));
	WIN_ASSERT_EQUAL(exact.NumMismatches, expectedNumMismatches, _T("Failed basespace conversion of the exact NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCAATGAGCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	expectedNumMismatches = 1;
	WIN_ASSERT_EQUAL(oneSNP.Reference, expectedReference, _T("Failed basespace conversion of the oneSNP reference.\n"));
	WIN_ASSERT_EQUAL(oneSNP.Query, expectedQuery, _T("Failed basespace conversion of the oneSNP query.\n"));
	WIN_ASSERT_EQUAL(oneSNP.BaseQualities, expectedQualities, _T("Failed basespace conversion of the oneSNP base qualities.\n"));
	WIN_ASSERT_EQUAL(oneSNP.NumMismatches, expectedNumMismatches, _T("Failed basespace conversion of the oneSNP NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	expectedNumMismatches = 0;
	WIN_ASSERT_EQUAL(oneError.Reference, expectedReference, _T("Failed basespace conversion of the oneError reference.\n"));
	WIN_ASSERT_EQUAL(oneError.Query, expectedQuery, _T("Failed basespace conversion of the oneError query.\n"));
	WIN_ASSERT_EQUAL(oneError.BaseQualities, expectedQualities, _T("Failed basespace conversion of the oneError base qualities.\n"));
	WIN_ASSERT_EQUAL(oneError.NumMismatches, expectedNumMismatches, _T("Failed basespace conversion of the oneError NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCAATGAGCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	expectedNumMismatches = 1;
	WIN_ASSERT_EQUAL(snpAndError.Reference, expectedReference, _T("Failed basespace conversion of the snpAndError reference.\n"));
	WIN_ASSERT_EQUAL(snpAndError.Query, expectedQuery, _T("Failed basespace conversion of the snpAndError query.\n"));
	WIN_ASSERT_EQUAL(snpAndError.BaseQualities, expectedQualities, _T("Failed basespace conversion of the snpAndError base qualities.\n"));
	WIN_ASSERT_EQUAL(snpAndError.NumMismatches, expectedNumMismatches, _T("Failed basespace conversion of the snpAndError NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCTTTGAGCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZabc";
	expectedNumMismatches = 2;
	WIN_ASSERT_EQUAL(errorAndSnp.Reference, expectedReference, _T("Failed basespace conversion of the errorAndSnp reference.\n"));
	WIN_ASSERT_EQUAL(errorAndSnp.Query, expectedQuery, _T("Failed basespace conversion of the errorAndSnp query.\n"));
	WIN_ASSERT_EQUAL(errorAndSnp.BaseQualities, expectedQualities, _T("Failed basespace conversion of the errorAndSnp base qualities.\n"));
	WIN_ASSERT_EQUAL(errorAndSnp.NumMismatches, expectedNumMismatches, _T("Failed basespace conversion of the errorAndSnp NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAG-TGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCAGCTGAGCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZabcd";
	expectedNumMismatches = 1;
	WIN_ASSERT_EQUAL(shortInsertion.Reference, expectedReference, _T("Failed basespace conversion of the shortInsertion reference.\n"));
	WIN_ASSERT_EQUAL(shortInsertion.Query, expectedQuery, _T("Failed basespace conversion of the shortInsertion query.\n"));
	WIN_ASSERT_EQUAL(shortInsertion.BaseQualities, expectedQualities, _T("Failed basespace conversion of the shortInsertion base qualities.\n"));
	WIN_ASSERT_EQUAL(shortInsertion.NumMismatches, expectedNumMismatches, _T("Failed basespace conversion of the shortInsertion NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAG----TGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCAGTACGTGAGCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg";
	expectedNumMismatches = 4;
	WIN_ASSERT_EQUAL(longInsertion.Reference, expectedReference, _T("Failed basespace conversion of the longInsertion reference.\n"));
	WIN_ASSERT_EQUAL(longInsertion.Query, expectedQuery, _T("Failed basespace conversion of the longInsertion query.\n"));
	WIN_ASSERT_EQUAL(longInsertion.BaseQualities, expectedQualities, _T("Failed basespace conversion of the longInsertion base qualities.\n"));
	WIN_ASSERT_EQUAL(longInsertion.NumMismatches, expectedNumMismatches, _T("Failed basespace conversion of the longInsertion NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCA-TGAGCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZab";
	expectedNumMismatches = 1;
	WIN_ASSERT_EQUAL(expectedReference,     shortDeletion.Reference,     _T("Failed basespace conversion of the shortDeletion reference.\n"));
	WIN_ASSERT_EQUAL(expectedQuery,         shortDeletion.Query,         _T("Failed basespace conversion of the shortDeletion query.\n"));
	WIN_ASSERT_EQUAL(expectedQualities,     shortDeletion.BaseQualities, _T("Failed basespace conversion of the shortDeletion base qualities.\n"));
	WIN_ASSERT_EQUAL(expectedNumMismatches, shortDeletion.NumMismatches, _T("Failed basespace conversion of the shortDeletion NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCA-TGAGCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXYZab";
	expectedNumMismatches = 1;
	WIN_ASSERT_EQUAL(expectedReference,     shortDeletionAndError.Reference,     _T("Failed basespace conversion of the shortDeletionAndError reference.\n"));
	WIN_ASSERT_EQUAL(expectedQuery,         shortDeletionAndError.Query,         _T("Failed basespace conversion of the shortDeletionAndError query.\n"));
	WIN_ASSERT_EQUAL(expectedQualities,     shortDeletionAndError.BaseQualities, _T("Failed basespace conversion of the shortDeletionAndError base qualities.\n"));
	WIN_ASSERT_EQUAL(expectedNumMismatches, shortDeletionAndError.NumMismatches, _T("Failed basespace conversion of the shortDeletionAndError NumMismatches.\n"));

	expectedReference = "GAGAGGCAGAGGTTGCAGTGAGCTGAGATT";
	expectedQuery     = "GAGAGGCAGAGGTTGCA----GCTGAGATT";
	expectedQualities = "AABCDEFGHIJKLMNOPQRSTUVWXY";
	expectedNumMismatches = 4;
	WIN_ASSERT_EQUAL(expectedReference,     longDeletion.Reference,     _T("Failed basespace conversion of the longDeletion reference.\n"));
	WIN_ASSERT_EQUAL(expectedQuery,         longDeletion.Query,         _T("Failed basespace conversion of the longDeletion query.\n"));
	WIN_ASSERT_EQUAL(expectedQualities,     longDeletion.BaseQualities, _T("Failed basespace conversion of the longDeletion base qualities.\n"));
	WIN_ASSERT_EQUAL(expectedNumMismatches, longDeletion.NumMismatches, _T("Failed basespace conversion of the longDeletion NumMismatches.\n"));
}
END_TEST

BEGIN_TEST(CColorspaceUtilities_ConvertReadBasespaceToPseudoColorspace) {
	CColorspaceUtilities csu;
	CMosaikString s              = "AATCTGAGAGGCAGAGGTTGCAGTGAGCTGAGATTGTGCC";
	const CMosaikString expected = "ATGGCGGGGATCGGGACACTCGCCGGTGCGGGTACCCTA";
	csu.ConvertReadBasespaceToPseudoColorspace(s);
	WIN_ASSERT_EQUAL(expected, s, _T("Failed the colorspace to pseudo-colorspace test.\n"));
}
END_TEST

BEGIN_TEST(CColorspaceUtilities_ConvertReadColorspaceToPseudoColorspace) {
	CColorspaceUtilities csu;
	CMosaikString s              = "03030300033332333-1111";
	const CMosaikString expected = "ATATATAAATTTTGTTTNCCCC";
	csu.ConvertReadColorspaceToPseudoColorspace(s);
	WIN_ASSERT_EQUAL(expected, s, _T("Failed the colorspace to pseudo-colorspace test.\n"));
}
END_TEST
