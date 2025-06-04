o2_add_executable(cluster-size
                    SOURCES determineClusterSize.cxx
                    PUBLIC_LINK_LIBRARIES O2::CommonUtils O2::CommonConstants O2::CommonDataFormat O2::CCDB O2::DataFormatsMID O2::MIDBase O2::MIDClustering O2::DataFormatsMCH
                    COMPONENT_NAME mid)