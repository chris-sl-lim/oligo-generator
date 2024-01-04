from flask import Flask, render_template, request, jsonify
from flask_socketio import SocketIO
import oligo_generator.models.generator as ogu

app = Flask(__name__)
socketio = SocketIO(app, cors_allowed_origins='*')
status = None


@app.route('/')
def home():
    return render_template('index.html')


@socketio.on('connect')
def handle_connect():
    print('Client connected')


@socketio.on('disconnect')
def handle_disconnect():
    print('Client disconnected')


@socketio.on('start_task')
@app.route('/generateSequences', methods=['POST'])
def create_generated_sequences():
    global status
    try:
        # Get the data from the web form
        data = request.get_json()

        data = parse_json_input(data)

        # Separate into variables
        #base_sequence = data['inputSequence']
        #num_changes = int(data['numChanges'])

        # Prepare oligo generator
        og = prepare_oligo_generator(data)

        # Run the oligo generator
        result = generate_sequences(og)

        # Jsonify results
        response = jsonify(result)
        
        # Return the result
        return response
    
    except Exception as e:
        return jsonify({"error": str(e)})
    

@app.route('/checkInputs', methods=['POST'])
def check_inputs():
    try:
        # Get the data from the web form
        data = request.get_json()

        data = parse_json_input(data)

        # Separate into variables
        #base_sequence = data['inputSequence']
        #num_changes = int(data['numChanges'])

        # Prepare oligo generator
        og = prepare_oligo_generator(data)

        # Pack results
        result = {
            "base_amino_acid_sequence": og.base_aa_seq,
            "fixed_positions": og.change_nt_vector
        }

        # Jsonify results
        response = jsonify(result)
        
        # Return the result
        return response
    
    except Exception as e:
        return jsonify({"error": str(e)})
    

def prepare_oligo_generator(data):

    # Unpack dictionary
    base_sequence = data['inputSequence']
    num_changes = data['numChanges']
    fixed_aa_positions = data['aaFixedPositions']
    fixed_nt_positions = data['ntFixedPositions']
    fixed_aa_codes = data['aaFixedCodes']
    restricted_aa_sequences = data['aaRestrictedSequences']
    restricted_nt_sequences = data['ntRestrictedSequences']

    # Create the generator
    o = ogu.oligo_generator(base_nt_seq=base_sequence)

    # Set the number of changes
    o.num_changes = num_changes

    # Set fixed amino acid positions
    for value in fixed_aa_positions:
        o.set_aa_pos(value, False)

    # Set fixed nucleotide positions
    for value in fixed_nt_positions:
        o.set_nt_pos(value, False)

    # Set fixed amino acid codes
    for value in fixed_aa_codes:
        o.set_fixed_aa(value)

    # Set restricted amino acid sequences
    o.restricted_aa_sequences = restricted_aa_sequences

    # Set restricted nucleotide sequences
    o.restriction_sites = restricted_nt_sequences

    # Return oligo generator
    return o
    

def generate_sequences(o):   

    # Create new amino acid sequences
    status = "Generating amino acid sequences"
    o.generate_aa_sequences(s_io = socketio)
    # socketio.emit('update_progress', {'progress': 50})

    # Create nucleotide sequences
    status = "Generating nucleotide sequences"
    o.generate_nt_sequences(s_io = socketio)
    # socketio.emit('update_progress', {'progress': 100})

    # Remove empty nucleotide sequences
    generated_nt_sequences = list(o.generated_nt_seq)
    non_empty_indices = [i for i,x in enumerate(generated_nt_sequences) if x]
    generated_aa_sequences = [list(o.generated_aa_seq)[idx] for idx in non_empty_indices]
    changes = [o.generated_nt_seq_change_attempts[idx] for idx in non_empty_indices]

    # Return the generated sequences
    return {"base_amino_acid_sequence": o.base_aa_seq,
            "nucleotide_sequences": generated_nt_sequences, 
            "amino_acid_sequences": generated_aa_sequences, 
            "changes": changes}


def parse_json_input(data):
        
    # Convert NumChanges to an integer
    data['numChanges'] = int(data['numChanges'])

    # Amino acid fixed positions (turn into integer list)
    if data['aaFixedPositions']:
        data['aaFixedPositions'] = [int(value) for value in data['aaFixedPositions'].split(',')]
    else:
        data['aaFixedPositions'] = []

    # Nucleotide fixed positions (turn into integer list)
    if data['ntFixedPositions']:
        data['ntFixedPositions'] = [int(value) for value in data['ntFixedPositions'].split(',')]
    else:
        data['ntFixedPositions'] = []

    # Amino acid fixed codes (turn into list)
    if data['aaFixedCodes']:
        data['aaFixedCodes'] = data['aaFixedCodes'].split(',')
    else:
        data['aaFixedCodes'] = []

    # Parse amino acid restricted sequences
    if data['aaRestrictedSequences']:
        data['aaRestrictedSequences'] = data['aaRestrictedSequences'].split('\n')
    else:
        data['aaRestrictedSequences'] = []     

    return data
    

if __name__ == '__main__':
    # app.run(debug=True)
    socketio.run(app, port=8000, host='0.0.0.0', debug=True)
