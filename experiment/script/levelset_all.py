import levelset_mnist
import levelset_perf
import levelset_perf3d
import vertebra64

def latex_float(f):
    if f == "--":
        return f
    float_str = "{0:.2e}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
    else:
        return float_str

def latex_dk(f):
    if f == "--":
        return f
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
    else:
        return float_str

times = dict()
times.update(levelset_mnist.run_all())
times.update(levelset_perf.run_all())
times.update(levelset_perf3d.run_all())
times.update(vertebra64.run_all())

def time_string(times, fc, name):
    str = "\n& {}".format(name) + \
    "\n& {}".format(latex_float(times.get('MNIST', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_float(times.get('Vert-64', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_float(times.get('S2D(0.01)', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_float(times.get('S2D(0.1)', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_float(times.get('S3D(0.01)', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_float(times.get('S3D(0.1)', dict()).get(fc, dict()).get(name, '--')))
    return str

def dk_string(times, fc, name):
    str = "\n& $d_K$" + \
    "\n& {}".format(latex_dk(times.get('MNIST', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_dk(times.get('Vert-64', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_dk(times.get('S2D(0.01)', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_dk(times.get('S2D(0.1)', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_dk(times.get('S3D(0.01)', dict()).get(fc, dict()).get(name, '--'))) + \
    "\n& {}".format(latex_dk(times.get('S3D(0.1)', dict()).get(fc, dict()).get(name, '--')))
    return str

# times = {"MNIST" : {"Freudenthal" : {"Ripser" : 0.015}}}
# print(times.get("Vert-64", dict()))

table_str = r"""
\begin{tabular}{|c|c||c|c|c|c|c|c|}
\hline
& & \multicolumn{1}{c|}{MNIST} & \multicolumn{1}{c|}{Vert-64} & \multicolumn{1}{c|}{S2D(0.01)} & \multicolumn{1}{c|}{S2D(0.1)} & \multicolumn{1}{c|}{S3D(0.01)} & \multicolumn{1}{c|}{S3D(0.1)}\\
\hline
\hline
\parbox[t]{2mm}{\multirow{7}{*}{\rotatebox[origin=c]{90}{Freudenthal}}}""" + \
dk_string(times, 'Freudenthal', 'dK') + "\n\\\\\cline{2-8}" + \
time_string(times, 'Freudenthal', 'Ripser') + "\n\\\\\cline{2-8}" + \
time_string(times, 'Freudenthal', 'Dionysus') + "\n\\\\\cline{2-8}" +  \
time_string(times, 'Freudenthal', 'Gudhi') + "\n\\\\\cline{2-8}" +\
time_string(times, 'Freudenthal', 'BATS(c)') + "\n\\\\\cline{2-8}" +\
time_string(times, 'Freudenthal', 'BATS(u,s)') + "\n\\\\\cline{2-8}" +\
time_string(times, 'Freudenthal', 'BATS(u,c)') + "\n" + r"\\\hline\hline" +\
"\n" + r"\parbox[t]{2mm}{\multirow{5}{*}{\rotatebox[origin=c]{90}{Cubical}}}" + \
dk_string(times, 'Cubical', 'dK') + "\n\\\\\cline{2-8}" + \
time_string(times, 'Cubical', 'Gudhi') + "\n\\\\\cline{2-8}" +\
time_string(times, 'Cubical', 'BATS(c)') + "\n\\\\\cline{2-8}" +\
time_string(times, 'Cubical', 'BATS(u,s)') + "\n\\\\\cline{2-8}" +\
time_string(times, 'Cubical', 'BATS(u,c)') + "\n" + r"\\\hline" + \
"\n" + r"\end{tabular}"


print(table_str)
