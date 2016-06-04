from flask import Flask, request, render_template, send_file, make_response
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from krig import *
import os
from bs4 import BeautifulSoup

app = Flask(__name__)


def kkk():
    # import and munge data
    APP_ROOT = os.path.dirname(os.path.abspath(__file__))
    folder = os.path.join(APP_ROOT, os.path.join('static', 'data'))
    z = open(os.path.join(folder, 'ZoneA.dat'), 'r').readlines()
    z = [i.strip().split() for i in z[10:]]
    z = np.array(z, dtype=np.float)
    z = DataFrame(z, columns=['x', 'y', 'thk', 'por', 'perm',
                              'lperm', 'lpermp', 'lpermr'])

    P = np.array(z[['x', 'y', 'por']])
    bw = 500  # bandwidth, plus or minus 250 meters
    # lags in 500 meter increments from zero to 10,000
    hs = np.arange(0, 10500, bw)
    sv = SV(P, hs, bw)
    vfunc_model = np.vectorize(spherical)
    sp = cvmodel(P, model=vfunc_model, hs=np.arange(0, 10500, 500), bw=500)

    X0, X1 = P[:, 0].min(), P[:, 0].max()
    Y0, Y1 = P[:, 1].min(), P[:, 1].max()
    Z = np.zeros((80, 100))
    V = np.zeros((80, 100))
    dx, dy = (X1 - X0)/100.0, (Y1 - Y0)/80.0

    model = cvmodel(P, vfunc_model, hs, bw)

    start = time.time()
    m = 80
    for i in xrange(m):
        percent = int(i / m * 100)
        print '[' + percent * '=' + (100-percent) * ' ' + '] ' +\
            str(percent) + '%'
        for j in xrange(100):
            Z[i, j], V[i, j] = krige(P, model, hs, bw, (dy*j, dx*i), 16)
    print '[' + 100 * '=' + '] ' + '100%' + '\nCompleted!'
    #        Z[i, j]=answer[0]
    #        V[i, j]=answer[1]
    print "\nTotal time: ", (time.time() - start)

    # Plot
    matplotlib.use('Agg')
    cdict = {'red':   ((0.0, 1.0, 1.0),
                       (0.5, 225/255., 225/255.),
                       (0.75, 0.141, 0.141),
                       (1.0, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),
                       (0.5, 57/255., 57/255.),
                       (0.75, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue':  ((0.0, 0.376, 0.376),
                       (0.5, 198/255., 198/255.),
                       (0.75, 1.0, 1.0),
                       (1.0, 0.0, 0.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(
                                                'my_colormap', cdict, 256)

    fig, ax = plt.subplots()
    H = np.zeros_like(Z)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            H[i, j] = np.round(Z[i, j] * 3)
    ax.matshow(H, cmap=my_cmap, interpolation='nearest')
    ax.scatter(z.x/200.0, z.y/200.0, facecolor='none', linewidths=0.75, s=50)
    ax.set_xlim(0, 99)
    ax.set_ylim(0, 80)
    ax.set_xticks([25, 50, 75], [5000, 10000, 15000])
    ax.set_yticks([25, 50, 75], [5000, 10000, 15000])
    fig.savefig('static/images/krigingpurple.svg', fmt='svg', dpi=200)

    return 'static/images/krigingpurple.svg'


@app.route('/')
def index():
    user_agent = request.headers.get('User-Agent')
    return render_template('index.html', outlook="Yuhao")


@app.route('/cal')
def cal():
    img_str = kkk()
    # return send_static_file(img_str, mimetype='image/svg+xml')
    # res = make_response(img_str)
    # res.headers["content-type"] = "text/plain"
    # return res
    data = ""
    with open(img_str, 'r') as file:
        for line in file.readlines():
            data += line.rstrip('\n')
    soup = BeautifulSoup(data, "html.parser")
    returnData = str(soup.svg)
    res = make_response(returnData)
    res.headers["content-type"] = "text/plain"
    return res

if __name__ == '__main__':
    app.run(debug=True, port=5000)
