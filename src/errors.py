class AuthenticationError(Exception):
    pass

class UserAlreadyExistsError(Exception):
    pass

class ValidationError(Exception):
    pass

class DatabaseConnectionError(Exception):
    pass