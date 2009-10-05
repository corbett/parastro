#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "vtkTestDriver.h"



/* Forward declare test functions. */
int otherCoordinate(int, char*[]);
int TestPriorityStreaming(int, char*[]);
int LoadOpenGLExtension(int, char*[]);
int SurfacePlusEdges(int, char*[]);
int TestActorLightingFlag(int, char*[]);
int TestAnimationScene(int, char*[]);
int TestBlurAndSobelPasses(int, char*[]);
int TestDynamic2DLabelMapper(int, char*[]);
int TestFBO(int, char*[]);
int TestGaussianBlurPass(int, char*[]);
int TestGradientBackground(int, char*[]);
int TestInteractorTimers(int, char*[]);
int TestLabelPlacer(int, char*[]);
int TestLabelPlacer2D(int, char*[]);
int TestLabelPlacerCoincidentPoints(int, char*[]);
int TestLabelPlacementMapper(int, char*[]);
int TestLabelPlacementMapper2D(int, char*[]);
int TestLabelPlacementMapperCoincidentPoints(int, char*[]);
int TestLightActor(int, char*[]);
int TestLODActor(int, char*[]);
int TestManyActors(int, char*[]);
int TestOrderedTriangulator(int, char*[]);
int TestOpacity(int, char*[]);
int TestOSConeCxx(int, char*[]);
int TestPOVExporter(int, char*[]);
int TestResetCameraVerticalAspectRatio(int, char*[]);
int TestResetCameraVerticalAspectRatioParallel(int, char*[]);
int TestSobelGradientMagnitudePass(int, char*[]);
int TestShadowMapPass(int, char*[]);
int TestTextActorAlphaBlending(int, char*[]);
int TestTextActorDepthPeeling(int, char*[]);
int TestTextActor3DAlphaBlending(int, char*[]);
int TestTextActor3DDepthPeeling(int, char*[]);
int TestTilingCxx(int, char*[]);
int TestTranslucentLUTAlphaBlending(int, char*[]);
int TestTranslucentLUTDepthPeeling(int, char*[]);
int TestTranslucentLUTDepthPeelingPass(int, char*[]);
int TestTranslucentLUTTextureAlphaBlending(int, char*[]);
int TestTranslucentLUTTextureDepthPeeling(int, char*[]);
int TestGenericVertexAttributesGLSLCxx(int, char*[]);
int TestGenericVertexAttributesGLSLAlphaBlending(int, char*[]);
int TestGenericVertexAttributesGLSLDepthPeelingPass(int, char*[]);


/* Create map.  */

typedef int (*MainFuncPointer)(int , char*[]);
typedef struct
{
  const char* name;
  MainFuncPointer func;
} functionMapEntry;

functionMapEntry cmakeGeneratedFunctionMapEntries[] = {
    {
    "otherCoordinate",
    otherCoordinate
  },
  {
    "TestPriorityStreaming",
    TestPriorityStreaming
  },
  {
    "LoadOpenGLExtension",
    LoadOpenGLExtension
  },
  {
    "SurfacePlusEdges",
    SurfacePlusEdges
  },
  {
    "TestActorLightingFlag",
    TestActorLightingFlag
  },
  {
    "TestAnimationScene",
    TestAnimationScene
  },
  {
    "TestBlurAndSobelPasses",
    TestBlurAndSobelPasses
  },
  {
    "TestDynamic2DLabelMapper",
    TestDynamic2DLabelMapper
  },
  {
    "TestFBO",
    TestFBO
  },
  {
    "TestGaussianBlurPass",
    TestGaussianBlurPass
  },
  {
    "TestGradientBackground",
    TestGradientBackground
  },
  {
    "TestInteractorTimers",
    TestInteractorTimers
  },
  {
    "TestLabelPlacer",
    TestLabelPlacer
  },
  {
    "TestLabelPlacer2D",
    TestLabelPlacer2D
  },
  {
    "TestLabelPlacerCoincidentPoints",
    TestLabelPlacerCoincidentPoints
  },
  {
    "TestLabelPlacementMapper",
    TestLabelPlacementMapper
  },
  {
    "TestLabelPlacementMapper2D",
    TestLabelPlacementMapper2D
  },
  {
    "TestLabelPlacementMapperCoincidentPoints",
    TestLabelPlacementMapperCoincidentPoints
  },
  {
    "TestLightActor",
    TestLightActor
  },
  {
    "TestLODActor",
    TestLODActor
  },
  {
    "TestManyActors",
    TestManyActors
  },
  {
    "TestOrderedTriangulator",
    TestOrderedTriangulator
  },
  {
    "TestOpacity",
    TestOpacity
  },
  {
    "TestOSConeCxx",
    TestOSConeCxx
  },
  {
    "TestPOVExporter",
    TestPOVExporter
  },
  {
    "TestResetCameraVerticalAspectRatio",
    TestResetCameraVerticalAspectRatio
  },
  {
    "TestResetCameraVerticalAspectRatioParallel",
    TestResetCameraVerticalAspectRatioParallel
  },
  {
    "TestSobelGradientMagnitudePass",
    TestSobelGradientMagnitudePass
  },
  {
    "TestShadowMapPass",
    TestShadowMapPass
  },
  {
    "TestTextActorAlphaBlending",
    TestTextActorAlphaBlending
  },
  {
    "TestTextActorDepthPeeling",
    TestTextActorDepthPeeling
  },
  {
    "TestTextActor3DAlphaBlending",
    TestTextActor3DAlphaBlending
  },
  {
    "TestTextActor3DDepthPeeling",
    TestTextActor3DDepthPeeling
  },
  {
    "TestTilingCxx",
    TestTilingCxx
  },
  {
    "TestTranslucentLUTAlphaBlending",
    TestTranslucentLUTAlphaBlending
  },
  {
    "TestTranslucentLUTDepthPeeling",
    TestTranslucentLUTDepthPeeling
  },
  {
    "TestTranslucentLUTDepthPeelingPass",
    TestTranslucentLUTDepthPeelingPass
  },
  {
    "TestTranslucentLUTTextureAlphaBlending",
    TestTranslucentLUTTextureAlphaBlending
  },
  {
    "TestTranslucentLUTTextureDepthPeeling",
    TestTranslucentLUTTextureDepthPeeling
  },
  {
    "TestGenericVertexAttributesGLSLCxx",
    TestGenericVertexAttributesGLSLCxx
  },
  {
    "TestGenericVertexAttributesGLSLAlphaBlending",
    TestGenericVertexAttributesGLSLAlphaBlending
  },
  {
    "TestGenericVertexAttributesGLSLDepthPeelingPass",
    TestGenericVertexAttributesGLSLDepthPeelingPass
  },

  {0,0}
};

/* Allocate and create a lowercased copy of string
   (note that it has to be free'd manually) */

char* lowercase(const char *string)
{
  char *new_string, *p;

#ifdef __cplusplus
  new_string = static_cast<char *>(malloc(sizeof(char) *
    static_cast<size_t>(strlen(string) + 1)));
#else
  new_string = (char *)(malloc(sizeof(char) * (size_t)(strlen(string) + 1)));
#endif

  if (!new_string)
    {
    return 0;
    }
  strcpy(new_string, string);
  p = new_string;
  while (*p != 0)
    {
#ifdef __cplusplus
    *p = static_cast<char>(tolower(*p));
#else
    *p = (char)(tolower(*p));
#endif

    ++p;
    }
  return new_string;
}

int main(int ac, char *av[])
{
  int i, NumTests, testNum, partial_match;
  char *arg, *test_name;
  int count;
  int testToRun = -1;

  
    
  for(count =0; cmakeGeneratedFunctionMapEntries[count].name != 0; count++)
    {
    }
  NumTests = count;
  /* If no test name was given */
  /* process command line with user function.  */
  if (ac < 2)
    {
    /* Ask for a test.  */
    printf("Available tests:\n");
    for (i =0; i < NumTests; ++i)
      {
      printf("%3d. %s\n", i, cmakeGeneratedFunctionMapEntries[i].name);
      }
    printf("To run a test, enter the test number: ");
    fflush(stdout);
    testNum = 0;
    if( scanf("%d", &testNum) != 1 )
      {
      printf("Couldn't parse that input as a number\n");
      return -1;
      }
    if (testNum >= NumTests)
      {
      printf("%3d is an invalid test number.\n", testNum);
      return -1;
      }
    testToRun = testNum;
    ac--;
    av++;
    }
  partial_match = 0;
  arg = 0;
  /* If partial match is requested.  */
  if(testToRun == -1 && ac > 1)
    {
    partial_match = (strcmp(av[1], "-R") == 0) ? 1 : 0;
    }
  if (partial_match && ac < 3)
    {
    printf("-R needs an additional parameter.\n");
    return -1;
    }
  if(testToRun == -1)
    {
    arg = lowercase(av[1 + partial_match]);
    }
  for (i =0; i < NumTests && testToRun == -1; ++i)
    {
    test_name = lowercase(cmakeGeneratedFunctionMapEntries[i].name);
    if (partial_match && strstr(test_name, arg) != NULL)
      {
      testToRun = i;
      ac -=2;
      av += 2;
      }
    else if (!partial_match && strcmp(test_name, arg) == 0)
      {
      testToRun = i;
      ac--;
      av++;
      }
    free(test_name);
    }
  if(arg)
    {
    free(arg);
    }
  if(testToRun != -1)
    {
    int result;
vtkFloatingPointExceptions::Enable();
    try {
    result = (*cmakeGeneratedFunctionMapEntries[testToRun].func)(ac, av);
    }
    catch(vtkstd::exception& e)
      {
      fprintf(stderr, "Test driver caught exception: [%s]\n", e.what());
      result = -1;
      }
    return result;
    }
  
  
  /* Nothing was run, display the test names.  */
  printf("Available tests:\n");
  for (i =0; i < NumTests; ++i)
    {
    printf("%3d. %s\n", i, cmakeGeneratedFunctionMapEntries[i].name);
    }
  printf("Failed: %s is an invalid test name.\n", av[1]);
  
  return -1;
}
