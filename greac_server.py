import os
import subprocess

from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/run-script', methods=['POST'])
def run_script():
    data = request.json

    train_dir = data.get('train_dir')
    test_dir = data.get('test_dir')
    group = data.get('group_name')
    window = data.get('window_size')
    metric = data.get('metric')
    cache = data.get('cache')

    print(data)
    script_path = f'./scripts/local/benchmark.sh'

    if not os.path.exists(script_path):
        return jsonify({"error": "Script não encontrado"}), 404

    try:
        result = subprocess.run(
            [
                script_path,
                train_dir,
                test_dir,
                group,
                window,
                metric,
                cache
            ],
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
