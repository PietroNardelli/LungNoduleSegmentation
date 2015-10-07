cmake_minimum_required(VERSION 2.8.9)

project(LungNoduleSegmentation)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/4.3/Extensions/LungNoduleSegmentation")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Pietro Nardelli (University College Cork)")
set(EXTENSION_DESCRIPTION "Scripted module for lung segmentation")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/LungNoduleSegmentation/LungNoduleSegmentation.png?revision=21747&view=co")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/e/e2/Slicer-r19441-LungNoduleSegmentation-screenshot.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(LungNoduleSegmentation)
add_subdirectory(LungNoduleSegmentationCLI)
add_subdirectory(LungSegmentationCLI)
add_subdirectory(LungFissureSegmentationCLI)
add_subdirectory(LungLobeSegmentationCLI)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})