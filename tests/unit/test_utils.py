from os import environ as env
from pathlib import Path
from vac2fost.vac2fost import shell_path

class TestShellPathReader:

    def test_home(self):
        assert shell_path("$HOME") == Path.home()

    def test_underscore(self):
        env["MCFOST_NULL"] = "notarealvar"
        assert shell_path("$HOME/$MCFOST_NULL") == Path.home()/"notarealvar"

    def test_double_underscore(self):
        env["MCFOST_extra_NULL"] = "notarealvar"
        assert shell_path("$HOME/$MCFOST_extra_NULL") == Path.home()/"notarealvar"

    def test_digits(self):
        env["MCFOST_002NotIntersing"] = "notarealvar"
        assert shell_path("$HOME/$MCFOST_002NotIntersing") == Path.home()/"notarealvar"
