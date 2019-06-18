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
