from flask import Flask, render_template, request, jsonify
from threading import Thread
import oligo_generator.models.generator as ogu

app = Flask(__name__)
status = None

@app.route('/')
def home():
    return render_template('index.html')

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

    # Return oligo generator
    return o
    

def generate_sequences(o):

    # Create new amino acid sequences
    status = "Generating amino acid sequences"
    o.generate_aa_sequences()

    # Create nucleotide sequences
    status = "Generating nucleotide sequences"
    o.generate_nt_sequences()

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

    # Turn the aaFixedPosition and ntFixedPositions into a list
    if data['aaFixedPositions']:
        data['aaFixedPositions'] = [int(value) for value in data['aaFixedPositions'].split(',')]
    else:
        data['aaFixedPositions'] = []

    if data['ntFixedPositions']:
        data['ntFixedPositions'] = [int(value) for value in data['ntFixedPositions'].split(',')]
    else:
        data['ntFixedPositions'] = []

    return data
    

if __name__ == '__main__':
    app.run(debug=True)
