project "Catch2"
	location "../../build"
	kind "StaticLib"
	language "C++"
	cppdialect "C++14"
	staticruntime "on"
	systemversion "latest"

	targetdir ("../../build/bin/" .. outputdir .. "%{prj.name}")
	objdir ("../../build/bin-int/" .. outputdir .. "%{prj.name}")

	files
	{
		"catch_amalgamated.hpp",
		"catch_amalgamated.cpp",
	}

	filter "configurations:Debug"
		runtime "Debug"
		symbols "on"

	filter "configurations:Release"
		runtime "Release"
		optimize "on"

	filter "configurations:Distribute"
		runtime "Release"
		optimize "on"