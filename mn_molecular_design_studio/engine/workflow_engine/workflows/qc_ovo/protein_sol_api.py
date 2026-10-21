PROTEIN_SOL_PATH = "/opt/protein-sol-sequence-prediction-software"
import os
import platform
import shutil
import subprocess
import tempfile

from typing import List


# def protein_sol(event, context):
#     sequences = event.get('sequences')
#     if not sequences:
#         return {
#             'statusCode': 400,
#             'body': json.dumps(
#                 {
#                     'status': 'Error'
#                 }
#             )
#         }

#     predictions = call_protein_sol(sequences)

#     return {
#         'statusCode': 200,
#         'body': json.dumps(
#             {'predictions': predictions}
#         )
#     }


def protein_sol(sequences: List[str]):
    if not sequences:
        raise ValueError("No sequences provided")

    predictions = call_protein_sol(sequences)
    return predictions


def call_protein_sol(sequences):
    with tempfile.TemporaryDirectory() as tmpdirname:
        copy_script_files(tmpdirname)
        process_sequences(sequences, tmpdirname)
        return parse_predictions_file(os.path.join(tmpdirname, "seq_prediction.txt"))


def copy_script_files(dest):
    protein_sol_path = PROTEIN_SOL_PATH
    script_filenames = [
        "multiple_prediction_wrapper_export.sh",
        "fasta_seq_reformat_export.pl",
        "seq_compositions_perc_pipeline_export.pl",
        "server_prediction_seq_export.pl",
        "seq_props_ALL_export.pl",
        "profiles_gather_export.pl",
        "ss_propensities.txt",
        "seq_reference_data.txt",
    ]

    for filename in script_filenames:
        shutil.copyfile(os.path.join(protein_sol_path, filename), os.path.join(dest, filename))

    os.chmod(os.path.join(dest, "multiple_prediction_wrapper_export.sh"), 0o775)


def process_sequences(sequences, directory):
    with tempfile.NamedTemporaryFile(mode="w+") as sequence_file:
        for i, sequence in enumerate(sequences):
            sequence_file.write(f">header{i}\n{sequence}\n")
        sequence_file.seek(0)
        if platform.system() == "Windows":
            subprocess.run(
                [
                    "C:\Program Files\Git\git-bash.exe",
                    os.path.join(directory, "multiple_prediction_wrapper_export.sh"),
                    sequence_file.name,
                ],
                cwd=directory,
            )
        else:
            subprocess.run(["./multiple_prediction_wrapper_export.sh", sequence_file.name], cwd=directory)


def parse_predictions_file(path):
    with open(path) as prediction_file:
        parser = SeqPredictionParser(prediction_file)
        predictions = parser.parse()
        return predictions


class SeqCompositionsParser:
    def __init__(self, file):
        self.lines = file.read().splitlines()
        self._current_line_number = 15
        self.headers = None

    def parse(self):
        headers = self.parse_headers()
        compositions = self.parse_compositions(headers)
        sums = self.parse_sums(self._current_line_number)
        return (compositions, sums)

    def parse_headers(self):
        headers_strings = (self.lines[4], self.lines[6], self.lines[8], self.lines[10], self.lines[12])
        self.headers = [[x.strip() for x in header_line] for header_line in [x.split(",") for x in headers_strings]]
        return self.headers

    def parse_compositions(self, headers):
        res = []
        while not self.lines[self._current_line_number].startswith("SUMS"):
            data = [x.split(",") for x in self.lines[self._current_line_number : self._current_line_number + 5]]
            res_single = {}
            for j, data_line in enumerate(data):
                data_line_stripped = [x.strip() for x in data_line]
                res_single[headers[j][0]] = dict(zip(headers[j][1:], data_line_stripped[1:]))
            res.append(res_single)
            self._current_line_number += 6
        return res

    def parse_sums(self, line_number):
        sums = {}
        sums["K-R"], sums["D-E"] = self.lines[line_number].split(" = ")[1].split()
        sums["K"], sums["R"], sums["D"], sums["E"] = self.lines[line_number + 2].split(" = ")[1].split()
        return sums


class SeqPredictionParser:
    def __init__(self, file):
        self.lines = file.read().splitlines()
        self._current_line_number = 14

    def parse(self):
        headers = self.parse_headers()
        return self.parse_prediction(headers)

    def parse_headers(self):
        headers_strings = self.lines[10:13]
        self.headers = [x.split(",") for x in headers_strings]
        return self.headers

    def parse_prediction(self, headers):
        res = []

        while self._current_line_number < len(self.lines):
            data = [x.split(",") for x in self.lines[self._current_line_number : self._current_line_number + 3]]
            res_single = {}

            predictions = [x.strip() for x in data[0]]
            res_single.update(zip(headers[0][2:], [float(x) for x in predictions[2:]]))

            deviations = [x.strip() for x in data[2]]
            res_single["deviations"] = dict(zip(headers[2][2:], [float(x) for x in deviations[2:]]))

            self._current_line_number += 3

            for i in range(2, 4):
                line_to_process = self.lines[self._current_line_number + i].split()[-1].split(",")
                res_single[line_to_process[0]] = ([None] * 10) + [float(x) for x in line_to_process[2:]] + ([None] * 10)

            res.append(res_single)
            self._current_line_number += 5

        return res


if __name__ == "__main__":
    protein_sol(["ASDASD"])
