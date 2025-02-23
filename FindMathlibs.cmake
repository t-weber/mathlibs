#
# mathlibs
# @author Tobias Weber (orcid: 0000-0002-7230-1932)
# @date 22-feb-2025
# @license see 'LICENSE' file
#

find_path(Mathlibs_INCLUDE_DIR
	NAMES mathlibs/tensor.h mathlibs/tensor_stat.h mathlibs/variadic_algos.h
	HINTS /usr/local/include
	DOC "mathlibs"
)


find_package_handle_standard_args(Mathlibs
	FOUND_VAR Mathlibs_FOUND
	REQUIRED_VARS Mathlibs_INCLUDE_DIR
)


if(Mathlibs_FOUND)
	set(Mathlibs_INCLUDE_DIRS ${Mathlibs_INCLUDE_DIR})
	message("Mathlibs include directory: ${Mathlibs_INCLUDE_DIRS}.")
else()
	message("Mathlibs was not found.")
endif()
