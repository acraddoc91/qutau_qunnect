// qutau_dll.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"

#define DLLEXPORT extern "C" __declspec(dllexport)

#include <stdio.h>
#include "tdcbase.h"


DLLEXPORT Int32 TDC_generateTimestamps_flat(double * par,	Int32 count) {
	Int32 rc;
	rc = TDC_generateTimestamps(SIM_FLAT, par, count);
	return rc;
}

DLLEXPORT Int32 TDC_generateTimestamps_norm(double * par, Int32 count) {
	Int32 rc;
	rc = TDC_generateTimestamps(SIM_NORMAL, par, count);
	return rc;
}