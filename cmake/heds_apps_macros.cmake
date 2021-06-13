# create information for future target creation based on user provided metadata
macro(disable_usual_warnings _name)
  IF(WIN32)
    IF(NOT BORLAND)
    IF(NOT CYGWIN)
      IF(NOT MINGW)
        target_compile_definitions(${_name} PRIVATE
            _CRT_FAR_MAPPINGS_NO_DEPRECATE
            _CRT_IS_WCTYPE_NO_DEPRECATE
            _CRT_MANAGED_FP_NO_DEPRECATE
            _CRT_NONSTDC_NO_DEPRECATE
            _CRT_SECURE_NO_DEPRECATE
            _CRT_SECURE_NO_DEPRECATE_GLOBALS
            _CRT_setERRORMODE_BEEP_SLEEP_NO_DEPRECATE
            _CRT_TIME_FUNCTIONS_NO_DEPRECATE
            _CRT_VCCLRIT_NO_DEPRECATE
            _SCL_SECURE_NO_DEPRECATE
            )     
        set(_BUILD_FLAGS /wd4251 /wd4267 /wd4244)
        target_compile_options(${_name} PRIVATE ${_BUILD_FLAGS})
      ENDIF(NOT MINGW)
    ENDIF(NOT CYGWIN)
    ENDIF(NOT BORLAND)
  ENDIF(WIN32)
endmacro()

