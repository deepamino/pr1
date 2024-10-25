import os
from flask import Flask, render_template, request, redirect, url_for
from werkzeug.utils import secure_filename
from process.translate_factory import TranslatorFactory
from process.api_manager import ApiManager

app = Flask(__name__)

UPLOAD_FOLDER = './uploads/'
SAVE_FOLDER = './results/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['SAVE_FOLDER'] = SAVE_FOLDER

@app.route('/deepamino/search')
def index():
    return render_template('index.html')

@app.route('/deepamino/process_search', methods=['GET'])
def process_search():
    search_term = request.args.get('search_term')

    if search_term:
        api_manager = ApiManager()
        protein_info = api_manager.reques_info(search_term)
        protein_object = api_manager.get_sequences(search_term)
        
        if protein_info:
            return render_template('results.html', information=[protein_info], proteins=[protein_object])
        else:
            return f"No se encontró información para el ID: {search_term}"
    
    return "No se ingresó ningún término de búsqueda."

@app.route('/deepamino/process_most_search', methods=['GET'])
def process_most_search():
    search_id = request.args.get('search_term')
    api_manager = ApiManager()
    filename = search_id + '.fasta'
    
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    save_path = os.path.join(app.config['SAVE_FOLDER'], filename)

    translator = TranslatorFactory.initialize_collector('Biopython')

    try:
        protein_objects = translator.process_fasta(file_path, save_path)
        information_objects = []

        for protein in protein_objects:
            info = api_manager.reques_info(protein.id)
            information_objects.append(info)
    except Exception as e:
            return f"Error al procesar el archivo: {str(e)}"
        
    return render_template('results.html', proteins=protein_objects, information=information_objects)
            

@app.route('/deepamino/upload', methods=['POST'])
def upload_file():
    api_manager = ApiManager()
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
            information_objects = []

            for protein in protein_objects:
                info = api_manager.reques_info(protein.id)
                information_objects.append(info)
            
        except Exception as e:
            return f"Error al procesar el archivo: {str(e)}"
        
        return render_template('results.html', proteins=protein_objects, information=information_objects)

@app.route('/deepamino/results')
def results():
    protein_objects = request.args.getlist('proteins')
    information_objects = request.args.getlist('information')
    print(information_objects)
    return render_template('results.html', proteins=protein_objects, information=information_objects)

if __name__ == '__main__':
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER)
    app.run(debug=True)
