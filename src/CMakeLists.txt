include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES   	boxjenkins.c
					brent.c
					conjgrad.c
					conv.c
					ctsa.c
					dist.c
					emle.c
					erfunc.c
					filter.c
					hsfft.c
					initest.c
					lls.c
					lnsrchmp.c
					matrix.c
					neldermead.c
					newtonmin.c
					nls.c
					optimc.c
					pdist.c
					polyroot.c
					pred.c
					real.c
					regression.c
					secant.c
					spectrum.c
					stats.c
					talg.c
					errors.c
					unitroot.c
					wtmath.c
					wavefilt.c
					wavelib.c
                    )

set(HEADER_FILES    boxjenkins.h   	
					brent.h		
					conjgrad.h
					conv.h
					ctsa.h
					dist.h
					emle.h
					erfunc.h
					filter.h
					hsfft.h
					initest.h
					lls.h
					lnsrchmp.h
					matrix.h
					neldermead.h
					newtonmin.h
					nls.h
					optimc.h
					pdist.h
					polyroot.h
					pred.h
					real.h
					regression.h
					secant.h
					spectrum.h
					stats.h
					talg.h
					errors.h
					unitroot.h
					wtmath.h
					wavefilt.h
					wavelib.h
					waveletarima.h
                    )

add_library(ctsalib STATIC ${SOURCE_FILES} ${HEADER_FILES})
        
set_property(TARGET ctsalib PROPERTY FOLDER "lib")

target_include_directories(ctsalib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})