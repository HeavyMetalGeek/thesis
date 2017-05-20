#!/usr/bin/python3
"""
Name:           file_strings.py
Last Modified:  05/20/2017
Python Version: 3.4.3

Description:
    Module containing strings (or functions returning strings) used for creating
    OpenFOAM files which would normally be cumbersome to include within
    functional code.

Contents:
    func of_header:
    	Returns string for an OpenFOAM file header

Usage:
	of_header(loc, obj_name)
		loc:
			Name of the OpenFOAM folder containing the file
		obj_name:
			Name of the object described within the file 

Notes:
"""

__all__ = ["of_header"]

from string import Template

def of_header(loc, obj_name):
	hdr = Template(r"""
	/*--------------------------------*- C++ -*----------------------------------*\
	| =========                 |                                                 |
	| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
	|  \\    /   O peration     | Version:  3.0.x                                 |
	|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
	|    \\/     M anipulation  |                                                 |
	\*---------------------------------------------------------------------------*/
	FoamFile
	{
	    version     2.0;
	    format      ascii;
	    class       volVectorField;
	    location    "$loc";
	    object      $obj;
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	""")

	return hdr.safe_substitute(loc=loc, obj=obj_name)