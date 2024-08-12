workspace "ves"
	location "build"
	configurations {"Debug", "Release", "Distribute"}
	platforms{"x64"}
	startproject "tests"
	filter { "platforms:x64" }
	architecture "x86_64"

outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"


include "third_party/Catch2"
include "tests"
include "ves"
include "third_party/ptp/ptp"
project("ptp")
	location "build"
	targetdir ("build/bin/" .. outputdir .. "%{prj.name}")
	objdir ("build/bin-int/" .. outputdir .. "%{prj.name}")