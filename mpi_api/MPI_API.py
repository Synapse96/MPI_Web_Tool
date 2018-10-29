from flask import Flask, jsonify, request
from flask_cors import CORS
from werkzeug.utils import secure_filename
import matlab.engine
import math

app = Flask(__name__)
CORS(app)
eng = matlab.engine.start_matlab()


class InvalidUsage(Exception):
    status_code = 400

    def __init__(self, message, status_code=None, payload=None):
        Exception.__init__(self)
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        rv = dict(self.payload or ())
        rv['message'] = self.message
        return rv


@app.route('/LeftMPI', methods=['GET', 'POST'])
def calculate_lmpi():
    file = request.files['imageA']
    filename = secure_filename(file.filename)
    if ".raw" not in filename:
        raise InvalidUsage('Incorrect image format, must be .raw format', status_code=400)
    sfp = app.config['UPLOAD_FOLDER'] + '/' + filename
    print('saving file', sfp)
    file.save(sfp)
    print('saved',)
    cell_array = eng.calculate_mpi(sfp)
    res_dict = {'MPI': list(filter(lambda x: not math.isnan(x), cell_array[1:]))}
    print(res_dict)
    return jsonify(res_dict), 200


@app.route('/RightMPI', methods=['GET', 'POST'])
def calculate_rmpi():
    file = request.files['imageA']
    filename = secure_filename(file.filename)
    if ".raw" not in filename:
        raise InvalidUsage('Incorrect image format, must be .raw format', status_code=400)
    a_wave = app.config['UPLOAD_FOLDER'] + '/' + filename
    file.save(a_wave)
    print(file)
    file = request.files['imageB']
    filename = secure_filename(file.filename)
    if ".raw" not in filename:
        raise InvalidUsage('Incorrect image format, must be .raw format', status_code=400)
    b_wave = app.config['UPLOAD_FOLDER'] + '/' + filename
    file.save(b_wave)
    print(file)
    arr_a = convert(eng.calculate_a_interval(a_wave)[0])
    arr_b = convert(eng.calculate_b_interval(b_wave)[0])
    res_dict = {'a_Interval': [], 'b_Interval': [], 'MPI': []}
    if len(arr_a) > len(arr_b):
        arr_b += [0.0]*(len(arr_a) - len(arr_b))
    else:
        arr_a += [0.0]*(len(arr_b) - len(arr_a))

    print(arr_a,arr_b)
    res_dict['a_Interval'] = arr_a
    res_dict['b_Interval'] = arr_b
    res_dict['MPI'] += arr_subtract(arr_a, arr_b)
    print(res_dict)
    return jsonify(res_dict), 200


@app.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    response = jsonify(error.to_dict())
    response.status_code = error.status_code
    return response


def convert(arr):
    conv = arr._data.tolist()
    return conv


def arr_subtract(arr_a, arr_b):
    result = []
    n = len(arr_a)
    i = 0
    while i < n:
        a = arr_a[i]
        b = arr_b[i]
        if b != 0 or a != 0:
            val = eng.calc_RMPI(a, b)
            if val == -1:
                result += [0]
            else:
                result += [val]
        else:
            result += [0]
        i = i + 1
    return result


if __name__ == '__main__':
    eng.addpath(app.root_path + '/MPI/LMPI')
    eng.addpath(app.root_path + '/MPI/RMPI/a_wave')
    eng.addpath(app.root_path + '/MPI/RMPI/b_wave')
    app.config['UPLOAD_FOLDER'] = app.root_path + '/temp'
    app.run()
