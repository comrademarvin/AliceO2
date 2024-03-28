o2_add_executable(hvlv-ccdb
                    SOURCES HVLV.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsDCS O2::CCDB O2::MCHGlobalMapping O2::MCHMappingImpl4 O2::MCHStatus
                    COMPONENT_NAME mch)