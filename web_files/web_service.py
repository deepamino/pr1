import os
from flask import Flask, render_template, request, redirect, url_for
from werkzeug.utils import secure_filename
from process.translate_factory import TranslatorFactory

app = Flask(__name__)

# Configuración de la carpeta para subir archivos
UPLOAD_FOLDER = './uploads/'
SAVE_FOLDER = './results/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['SAVE_FOLDER'] = SAVE_FOLDER

@app.route('/deepamino/search')
def index():
    return render_template('index.html')

@app.route('/deepamino/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return "No se ha seleccionado ningún archivo."
    
    file = request.files['file']
    
    if file.filename == '':
        return "No se ha seleccionado ningún archivo."

    if file:
        filename = secure_filename(file.filename)
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        save_path = os.path.join(app.config['SAVE_FOLDER'], filename)
        file.save(file_path)
        
        translator = TranslatorFactory.initialize_collector('Biopython')

        try:
            protein_objects = translator.process_fasta(file_path, save_path)
        except Exception as e:
            return f"Error al procesar el archivo: {str(e)}"
        
        return render_template('results.html', proteins=protein_objects)

@app.route('/deepamino/results')
def results():
    protein_objects = request.args.getlist('proteins')
    return render_template('results.html', proteins=protein_objects)

if __name__ == '__main__':
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER)
    app.run(debug=True)
