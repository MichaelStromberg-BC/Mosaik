// ***************************************************************************
// RegexUtilitiesTest.cpp - provides unit tests for CRegexUtilities.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <string>
#include "MosaikString.h"
#include "RegexUtilities.h"
#include "WinUnit.h"

using namespace std;

BEGIN_TEST(CRegexUtilities_ConvertQualities) {

	string fastaQualities1 = "65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65";
	string fastaQualities2 = "65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65\n";
	string fastaQualities3 = "65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65\r\n";
	string fastaQualities4 = "65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 ";
	string fastaQualities5 = "65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 \n";
	string fastaQualities6 = "65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 65 \r\n";

	const CMosaikString expected = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
	CMosaikString s;

	CRegexUtilities::ConvertQualities(fastaQualities1, s);
	WIN_ASSERT_EQUAL(s, expected, _T("Failed conversion on a list of base qualities.\n"));

	CRegexUtilities::ConvertQualities(fastaQualities2, s);
	WIN_ASSERT_EQUAL(s, expected, _T("Failed conversion on a list of base qualities with CR.\n"));

	CRegexUtilities::ConvertQualities(fastaQualities3, s);
	WIN_ASSERT_EQUAL(s, expected, _T("Failed conversion on a list of base qualities with CRLF.\n"));

	CRegexUtilities::ConvertQualities(fastaQualities4, s);
	WIN_ASSERT_EQUAL(s, expected, _T("Failed conversion on a list of base qualities with trailing space.\n"));

	CRegexUtilities::ConvertQualities(fastaQualities5, s);
	WIN_ASSERT_EQUAL(s, expected, _T("Failed conversion on a list of base qualities with trailing space and CR.\n"));

	CRegexUtilities::ConvertQualities(fastaQualities6, s);
	WIN_ASSERT_EQUAL(s, expected, _T("Failed conversion on a list of base qualities with trailing space and CRLF.\n"));
}
END_TEST
