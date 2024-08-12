project "tests"
	location "../build"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++14"
	staticruntime "on"
	systemversion "latest"

	targetdir ("../build/bin/" .. outputdir .. "%{prj.name}")
	objdir ("../build/bin-int/" .. outputdir .. "%{prj.name}")

	IncludeDir = {}
	IncludeDir["Catch2"] = "../third_party/Catch2"
	IncludeDir["ves"] = "../ves/src"

	files
	{
		"src/**.h",
		"src/**.hpp",
		"src/**.cpp",
	}

	includedirs
	{
		"src",
		"%{IncludeDir.Catch2}",
		"../third_party/",
		"%{IncludeDir.ves}",
		
	}

	links {
		"Catch2",
		"ptp"
	}

	filter "configurations:Debug"
		runtime "Debug"
		symbols "on"

	filter "configurations:Release"
		runtime "Release"
		optimize "on"