import os
import subprocess
import tomllib

# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        # out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        out = _minimal_ext_cmd(['git', 'describe', '--tags', '--dirty'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

# Return version string from pyproject.toml
def project_version(toml="pyproject.toml"):
    try:
        with open(toml, "rb") as f:
            ver = tomllib.load(f)['project']['version']
    except OSError:
        ver = "Unknown"

    return ver

# Return version string
def get_version():
    ver = git_version()
    if not ver == "Unknown":
        return ver
    ver = project_version()
    if not ver == "Unknown":
        return 'v' + ver
    return ver

# ver, rev, git, dty, = (get_version().split('-') + (['']*4))[:4]
# print(ver, rev, git, dty)
