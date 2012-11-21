
# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#libs/bioinfo-libs'
commons_path = '#libs/common-libs'

env = Environment(tools = ['default', 'packaging'],
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp',
                  CPPPATH = ['#', '#src', '#include', '/usr/local/include', '/usr/include/libxml2', bioinfo_path, commons_path ],
                  LIBPATH = ['/usr/lib', '/usr/local/lib', '#libs', '#libs/common-libs/', commons_path ],
                  LIBS = ['argtable2', 'common', 'config', 'cprops', 'xml2', 'z'],
                  LINKFLAGS = ['-fopenmp'])
                  
if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3'

env['objects'] = []


##### Targets

# Compile dependencies
SConscript(['%s/bioformats/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path
            ], exports = ['env', 'debug'])

# Create binaries and copy them to 'bin' folder
progs = Program('hpg-fastq', source = [Glob('src/*.c')], exports = ['env', 'debug', 'commons_path', 'bioinfo_path' ])

env.Install('#bin', ['hpg-fastq.conf'])

# Create Debian package
#if 'debian' in COMMAND_LINE_TARGETS:
#    SConscript("deb/SConscript", exports = ['env'] )
