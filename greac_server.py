import os
import subprocess

from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/run-script', methods=['POST'])
def run_script():
    data = request.json

    input_file = data.get('input')
    group = data.get('group')
    window = data.get('window')

    script_path = './scripts/local/benchmark.sh'

    if not os.path.exists(script_path):
        return jsonify({"error": "Script não encontrado"}), 404

    try:
        result = subprocess.run(script_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        return jsonify({
            "stdout": result.stdout,
            "stderr": result.stderr,
            "status": "success"
        })
    except subprocess.CalledProcessError as e:
        return jsonify({
            "stdout": e.stdout,
            "stderr": e.stderr,
            "status": "error"
        }), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5005)
