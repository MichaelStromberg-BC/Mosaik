// ***************************************************************************
// SequenceUtilitiesTest.cpp - provides unit tests for CSequenceUtilities.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <string>
#include "SequenceUtilities.h"
#include "WinUnit.h"

using namespace std;

BEGIN_TEST(CSequenceUtilities_Chomp) {
	char expected[32], unix[32], dos[32];
	sprintf_s(expected, 32, "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	sprintf_s(unix,     32, "ABCDEFGHIJKLMNOPQRSTUVWXYZ\n");
	sprintf_s(dos,      32, "ABCDEFGHIJKLMNOPQRSTUVWXYZ\r\n");

	CSequenceUtilities::Chomp(unix);
	WIN_ASSERT_ZERO(strcmp(unix, expected), _T("Failed the UNIX chomp test.\n"));

	CSequenceUtilities::Chomp(dos);
	WIN_ASSERT_ZERO(strcmp(dos, expected), _T("Failed the DOS chomp test.\n"));
}
END_TEST

BEGIN_TEST(CSequenceUtilities_GetReverseComplement) {
	char expected[32], s[32];
	sprintf_s(expected, 32, "GGGGAAAAAAAATTTATATAT");
	sprintf_s(s,        32, "ATATATAAATTTTTTTTCCCC");

	CSequenceUtilities::GetReverseComplement(s, strlen(s));
	WIN_ASSERT_ZERO(strcmp(s, expected), _T("Failed the reverse complement test.\n"));
}
END_TEST

BEGIN_TEST(CSequenceUtilities_LowercaseSequence) {
	string s              = "OSCAR";
	const string expected = "oscar";
	CSequenceUtilities::LowercaseSequence(s);
	WIN_ASSERT_EQUAL(s, expected, _T("Failed the lowercase test.\n"));
}
END_TEST

BEGIN_TEST(CSequenceUtilities_ReverseSequence) {
	char expected[32], s[32];
	sprintf_s(expected, 32, "OSCAR");
	sprintf_s(s,        32, "RACSO");
	CSequenceUtilities::ReverseSequence(s, strlen(s));
	WIN_ASSERT_ZERO(strcmp(s, expected), _T("Failed the reverse test.\n"));
}
END_TEST

BEGIN_TEST(CSequenceUtilities_UppercaseSequence) {
	string s              = "oscar";
	const string expected = "OSCAR";
	CSequenceUtilities::UppercaseSequence(s);
	WIN_ASSERT_EQUAL(s, expected, _T("Failed the uppercase test.\n"));
}
END_TEST
