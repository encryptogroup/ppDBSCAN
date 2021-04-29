
add_executable(ppDBSCAN_distance ppDBSCAN_distance.cpp common/ppdbscan_coordination.cpp common/CircuitWrapper.cpp common/CircuitWrapper.h)
target_link_libraries(ppDBSCAN_distance ABY::aby ENCRYPTO_utils::encrypto_utils)

add_executable(ppDBSCAN_grouping ppDBSCAN_grouping.cpp common/ppdbscan_coordination.cpp common/CircuitWrapper.cpp common/CircuitWrapper.h)
target_link_libraries(ppDBSCAN_grouping ABY::aby ENCRYPTO_utils::encrypto_utils)

