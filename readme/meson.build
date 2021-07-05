project('tiny_md', 'c',
        default_options: [ 'c_std=gnu11', 'warning_level=2' ]
)

cc = meson.get_compiler('c')

m_dep = cc.find_library('m')

headless = executable('tiny_md',
    [ 'tiny_md.c', 'core.c' ],
    dependencies: [ m_dep ]
)

gl_dep = cc.find_library('GL', required: false)
glu_dep = cc.find_library('GLU', required: false)
glut_dep = cc.find_library('glut', required: false)

viz = executable('viz',
    [ 'viz.c', 'core.c' ],
    dependencies: [ m_dep, gl_dep, glu_dep, glut_dep ]
)
