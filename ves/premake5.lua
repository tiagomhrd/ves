project "ves"
	location "../build"
	kind "StaticLib"
	language "C++"
	cppdialect "C++14"
	staticruntime "on"

	targetdir ("../build/bin/" .. outputdir .. "%{prj.name}")
	objdir ("../build/bin-int/" .. outputdir .. "%{prj.name}")

	files
	{
		"src/**.h",
		"src/**.hpp",
		"src/**.cpp",
	}

	includedirs
	{
		"src/",
		"../third_party",
	}

	filter "configurations:Debug"
		runtime "Debug"
		symbols "on"
		defines {
			"VES_DEBUG"
		}

	filter "configurations:Release"
		runtime "Release"
		optimize "on"
		defines {
			"VES_RELEASE"
		}

	filter "configurations:Distribute"
		runtime "Release"
		optimize "on"
		defines {
			"VES_DISTRIBUTE"
		}
	