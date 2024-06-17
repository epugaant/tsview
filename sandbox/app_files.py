
from flask import request, jsonify
from flask import Flask,send_file,send_from_directory

app = Flask(__name__)
app.config["DEBUG"] = True

# The absolute path of the directory containing images for users to download
app.config["CLIENT_GAIA"] = "/Users/epuga/ESDC/TSViz/data/gaia/anonymous1690191210843_fits"

# The absolute path of the directory containing CSV files for users to download
app.config["CLIENT_JWST"] = "E:/AudiotoText/Flask_File_Downloads/filedownload/files/csv"


@app.route('/get-gaia/<gaia_filename>',methods = ['GET','POST'])
def get_csv(csv_filename):

    try:
        return send_from_directory(app.config["CLIENT_GAIA"], filename= gaia_filename, as_attachment=True)
    except FileNotFoundError:
        abort(404)


@app.route('/get-jwst/<jwst_filename>',methods = ['GET','POST'])
def get_pdf(pdf_filename):

    try:
        return send_from_directory(app.config["CLIENT_JWST"], filename=pdf_filename, as_attachment=True)
    except FileNotFoundError:
        abort(404)



def api_id():
    # Check if an ID was provided as part of the URL.
    # If ID is provided, assign it to a variable.
    # If no ID is provided, display an error in the browser.
    if 'id' in request.args:
        id = int(request.args['id'])
    else:
        return "Error: No id field provided. Please specify an id."

    # Create an empty list for our results
    results = []

    # Loop through the data and match results that fit the requested ID.
    # IDs are unique, but other fields might return many results
    for lc in lcs:
        if lc['id'] == id:
            results.append(lc)

    # Use the jsonify function from Flask to convert our list of
    # Python dictionaries to the JSON format.
    return jsonify(results)

if __name__ == '__main__':
    app.run(host='0.0.0.0',port = 8000, threaded = True, debug = True)

#http://0.0.0.0:8000/api/v1/resources/books?id=2
