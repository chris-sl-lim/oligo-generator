from flask import Flask, render_template, request, jsonify
import oligo_generator.models.generator as ogu

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/generateSequences', methods=['POST'])
def create_generated_sequences():
    try:
        # Get the data from the web form
        data = request.get_json()

        # Separate into variables
        base_sequence = data['inputSequence']
        num_changes = int(data['numChanges'])

        # Run the oligo generator
        result = generate_sequences(base_sequence, num_changes)
        
        # Return the result
        return jsonify(result)
    
    except Exception as e:
        return jsonify({"error": str(e)})
    

def generate_sequences(base_sequence, num_changes):

    # Create the generator
    o = ogu.oligo_generator(base_nt_seq=base_sequence)

    # Set the number of changes
    o.num_changes = num_changes

    # Create new amino acid sequences
    o.generate_aa_sequences()

    # Create nucleotide sequences
    o.generate_nt_sequences()

    # Return the generated sequences
    return {"base_amino_acid_sequence": o.base_aa_seq,
            "nucleotide_sequences": list(o.generated_nt_seq), 
            "amino_acid_sequences": list(o.generated_aa_seq), 
            "changes": o.generated_nt_seq_change_attempts}
    

if __name__ == '__main__':
    app.run(debug=True)