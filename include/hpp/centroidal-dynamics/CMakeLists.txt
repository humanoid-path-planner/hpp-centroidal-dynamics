# Declare Headers
SET(${PROJECT_NAME}_HEADERS
  config.hh
  util.hh
  logger.hh
  solver_LP_abstract.hh
  solver_LP_qpoases.hh
  solver_LP_clp.hh
  centroidal_dynamics.hh
  stop-watch.hh
  )

INSTALL(FILES
  ${${PROJECT_NAME}_HEADERS}
  DESTINATION include/hpp/centroidal-dynamics
  )
