$(document).on('keypress', function(e) {
    if(e.which == 13) { // Enter key
        if ($('#main_ui').find('#btn_login').is(':visible')) {
            $('#btn_login').click();
        } else if ($('#main_ui').find('#btn_register').is(':visible')) {
            $('#btn_register').click();
        } else if ($('#main_ui').find('#btn_reset_password_initiate').is(':visible')) {
            $('#btn_reset_password_initiate').click();
        } else if ($('#main_ui').find('#btn_reset_password_final').is(':visible')) {
            $('#btn_reset_password_final').click();
        }
    }
});