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

    # Return the generated sequences. The model keeps these lists aligned when
    # nucleotide candidates are filtered out.
    return {"base_amino_acid_sequence": o.base_aa_seq,
            "nucleotide_sequences": list(o.generated_nt_seq),
            "amino_acid_sequences": list(o.generated_aa_seq),
            "changes": list(o.generated_nt_seq_change_attempts)}


def _parse_position_list(value):
    if not value:
        return []

    return [int(item.strip()) for item in value.split(',')
            if item.strip()]


def _parse_sequence_list(value):
    if not value:
        return []

    if isinstance(value, list):
        values = value
    else:
        values = value.split('\n')

    return [item.strip().upper() for item in values if item.strip()]


def parse_json_input(data):
        
    # Convert NumChanges to an integer
    data['numChanges'] = int(data['numChanges'])
    data['inputSequence'] = data['inputSequence'].strip().upper()

    # Amino acid and nucleotide fixed positions are entered as 1-based values.
    data['aaFixedPositions'] = _parse_position_list(data['aaFixedPositions'])
    data['ntFixedPositions'] = _parse_position_list(data['ntFixedPositions'])

    # Amino acid fixed codes (turn into list)
    data['aaFixedCodes'] = _parse_sequence_list(data['aaFixedCodes'].replace(',', '\n'))

    # Parse restricted sequences
    data['aaRestrictedSequences'] = _parse_sequence_list(
        data['aaRestrictedSequences']
    )
    data['ntRestrictedSequences'] = _parse_sequence_list(
        data['ntRestrictedSequences']
    )

    return data
    

if __name__ == '__main__':
    # app.run(debug=True)
    socketio.run(app, port=8000, host='0.0.0.0', debug=True)
