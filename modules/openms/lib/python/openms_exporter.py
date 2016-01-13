from xml.sax import make_parser, ContentHandler
from optparse import OptionParser


def main():
    (options, _) = _parse_args()
    with open(options.output, "w") as out:
        parser = make_parser()
        handler = _get_handler(options)(out)
        parser.setContentHandler(handler)
        parser.parse(open(options.input, "r"))


def _get_handler(options):
    return handlers[options.type]


class OpenMsContentHandler(ContentHandler):

    def __record_values(self, keys, attrs):
        for key in keys:
            setattr(self, key, attrs.get(key, None))

    def _get_values(self, keys):
        return [getattr(self, key, "") for key in keys]

    def _set_attributes(self, name, attrs):
        for element_name, element_attributes in self.record_values.iteritems():
            if name == element_name:
                self.__record_values(element_attributes, attrs)

    def _write_line(self, line):
        self.output.write(line)
        self.output.write("\n")

    def startElement(self, name, attrs):
        self._set_attributes(name, attrs)

    def _handleElement(self, name):
        pass

    def endElement(self, name):
        self._handleElement(name)
        self._set_attributes(name, {})

    def _write_row(self, col_keys):
        row_values = self._get_values(col_keys)
        row = "\t".join(row_values)
        self._write_line(row)


class FeatureHullHandler(OpenMsContentHandler):
    record_values = {
        "feature": ["id"],
        "convexhull": ["nr"],
        "pt": ["x", "y"]
    }

    def __init__(self, output):
        self.output = output

    def _handleElement(self, name):
        if name == "pt":
            self._write_point()

    def _write_point(self):
        col_keys = ["id", "nr", "x", "y"]
        self._write_row(col_keys)


class PeptideHandler(OpenMsContentHandler):
    record_values = {
        "IdentificationRun": ["search_engine"],
        "PeptideIdentification": ["score_type", "significance_threshold", "MZ", "RT"],
        "PeptideHit": ["score", "sequence", "charge"],
    }

    def __init__(self, output):
        self.output = output

    def _handleElement(self, name):
        if name == "PeptideHit":
            self._write_peptide()

    def _write_peptide(self):
        col_keys = ["score", "sequence", "score_type", "charge", "MZ", "RT"]
        self._write_row(col_keys)


handlers = {
    "peptide": PeptideHandler,
    "feature_hull": FeatureHullHandler,
}


def _parse_args():
    parser = OptionParser()
    parser.add_option("--input", dest="input")
    parser.add_option("--output", dest="output")
    parser.add_option("--type", dest="type", choices=["peptide", "feature_hull"])
    return parser.parse_args()

if __name__ == "__main__":
    main()
