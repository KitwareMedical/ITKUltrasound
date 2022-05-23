set(test_input_dir ${${itk-module}_SOURCE_DIR}/test/Input)
set(test_baseline_dir ${${itk-module}_SOURCE_DIR}/test/Baseline)

list(FIND ITK_WRAP_IMAGE_DIMS 3 wrap_3_index)
if(wrap_3_index GREATER -1)
    itk_python_add_test(NAME PythonBModeImageFilterTestTiming
      TEST_DRIVER_ARGS
      --compare ${ITK_TEST_OUTPUT_DIR}/PythonBModeImageFilterTestTiming.mha
        DATA{${test_baseline_dir}/itkBModeImageFilterTestTimingBaseline.mha}
      COMMAND PythonBModeImageFilterTestTiming.py
        DATA{${test_input_dir}/uniform_phantom_8.9_MHz.mha}
        ${ITK_TEST_OUTPUT_DIR}/PythonBModeImageFilterTestTiming.mha
      )
endif()

itk_python_add_test(NAME PythonSpectra1DAveragingImageFilterTest
  TEST_DRIVER_ARGS
  --compare
    DATA{${test_input_dir}/spectra_pixel_image.mha}
    ${ITK_TEST_OUTPUT_DIR}/PythonSpectra1DAveragingImageFilterTest.mha
  COMMAND PythonSpectra1DAveragingImageFilterTest.py
    -o ${ITK_TEST_OUTPUT_DIR}/PythonSpectra1DAveragingImageFilterTest.mha
    -i DATA{${test_input_dir}/spectra_pixel_image.mha}
    -i DATA{${test_input_dir}/spectra_pixel_image.mha}
    -i DATA{${test_input_dir}/spectra_pixel_image.mha}
  )

itk_python_add_test(NAME PythonSpectra1DNormalizeImageFilterTest
  TEST_DRIVER_ARGS
  --compare
    DATA{${test_baseline_dir}/itkSpectra1DNormalizeImageFilterTestBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DNormalizeImageFilterTest.mha
  COMMAND PythonSpectra1DNormalizeImageFilterTest.py
    DATA{${test_input_dir}/spectra_pixel_image.mha}
    DATA{${test_input_dir}/spectra_pixel_image.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DNormalizeImageFilterTest.mha
  )
  
itk_python_add_test(NAME PythonCurvilinearArraySpecialCoordinatesImageTemplateMethodsTest
  COMMAND PythonCurvilinearArraySpecialCoordinatesImageTemplateMethodsTest.py
  )

itk_python_add_test(NAME PythonBackscatterImageFilterTest0
  TEST_DRIVER_ARGS
  --compareIntensityTolerance 0.01
  --compare
    DATA{${test_baseline_dir}/itkBackscatterImageFilterTest0Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/PythonBackscatterImageFilterTest0.mha
  --compare
    DATA{${test_baseline_dir}/itkBackscatterImageFilterTest1Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/PythonBackscatterImageFilterTest1.mha
  --compare
    DATA{${test_baseline_dir}/itkBackscatterImageFilterTest2Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/PythonBackscatterImageFilterTest2.mha
  COMMAND PythonBackscatterImageFilterTest.py
    -i DATA{${test_input_dir}/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    -a ${ITK_TEST_OUTPUT_DIR}/PythonBackscatterImageFilterTest0.mha
    -s ${ITK_TEST_OUTPUT_DIR}/PythonBackscatterImageFilterTest1.mha
    -n ${ITK_TEST_OUTPUT_DIR}/PythonBackscatterImageFilterTest2.mha
  )
