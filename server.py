from flask import Flask, request, redirect, url_for, render_template, jsonify
import pandas as pd
import subprocess
import threading
import os

app = Flask(__name__)

# Global variable to store progress
progress = 0

def get_filtered_df(gene_ids):
    """Retrieve filtered DataFrame based on gene_ids."""
    try:
        file_path = f'{gene_ids}_CRISPR_tgts.csv'
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            return df
        else:
            raise FileNotFoundError(f"No CSV file found for {gene_ids}")
    except Exception as e:
        print(f"An error occurred while loading the DataFrame: {str(e)}")
        return pd.DataFrame()

def run_crispr_tool(gene_ids):
    """Run the CRISPR web tool and update progress."""
    global progress
    try:
        process = subprocess.Popen(
            ['python3', 'crispr_webtooltest.py', gene_ids],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, universal_newlines=True
        )

        for line in process.stdout:
            line = line.strip()
            print(f"Subprocess Output: {line}")
            if "Progress" in line:
                progress = int(line.split()[1].replace('%', ''))
                print(f"Updated Progress: {progress}%")

        process.wait()
        progress = 100

        stderr_output = process.stderr.read()
        if stderr_output:
            print(f"Subprocess Error: {stderr_output}")

    except subprocess.CalledProcessError as e:
        print(f"CRISPR-Webtool Error: {e.stderr}")
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/submit', methods=['POST'])
def submit():
    gene_ids = request.form.get('gene_ids')
    print(f"Received Gene IDs: {gene_ids}")

    thread = threading.Thread(target=run_crispr_tool, args=(gene_ids,))
    thread.start()

    return redirect(url_for('progress_page', gene_ids=gene_ids))

@app.route('/progress/<gene_ids>')
def progress_page(gene_ids):
    return render_template('progress.html', gene_ids=gene_ids)

@app.route('/progress_status')
def progress_status():
    global progress
    return jsonify({'progress': progress})

@app.route('/dataframe/<gene_ids>')
def display_dataframe(gene_ids):
    try:
        filtered_df = get_filtered_df(gene_ids)
        filtered_html = filtered_df.to_html(classes='table table-striped', index=True)
        return render_template('dataframe.html', table=filtered_html, gene_ids=gene_ids, df=filtered_df.to_dict('records'))
    except FileNotFoundError as e:
        return f"File not found: {e}", 404

@app.route('/select_row', methods=['POST'])
def select_row():
    row_index = int(request.form.get('row_index'))
    gene_ids = request.form.get('gene_ids')
    filtered_df = get_filtered_df(gene_ids)
    
    if row_index in filtered_df.index:
        selected_target = filtered_df.loc[row_index]
        selected_csv_file = f'{gene_ids}_selected_row.csv'
        selected_target.to_frame().T.to_csv(selected_csv_file, index=False)
        
        command = ['python3', 'process_selected_row.py', selected_csv_file]
        result = subprocess.run(command, capture_output=True, text=True)
        
        processed_output_file = selected_csv_file.replace('.csv', '_processed_output.txt')
        if os.path.exists(processed_output_file):
            redirect_url = url_for('show_processed_txt', filename=os.path.basename(processed_output_file))
            return jsonify({'status': 'success', 'redirect_url': redirect_url})
        else:
            return jsonify({'status': 'error', 'message': 'Processed output file not found.'})
    else:
        return jsonify({'status': 'error', 'message': 'Invalid row index'})

@app.route('/show_processed_txt/<filename>')
def show_processed_txt(filename):
    try:
        file_path = os.path.join('.', filename)
        with open(file_path, 'r') as file:
            content = file.read()
        return render_template('show_txt.html', content=content)
    except FileNotFoundError:
        return "File not found", 404

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
