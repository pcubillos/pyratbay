from pygments.style import Style
from pygments.styles.friendly import FriendlyStyle
from pygments.token import Token


class CustomStyle(Style):
    styles = FriendlyStyle.styles.copy()

    background_color = "#f0f8ff"
    background_color = "#f0f2f6"
    styles[Token.Background] = background_color

    comment = "italic #bb3322"
    styles[Token.Comment] = comment
    styles[Token.Comment] = comment
    styles[Token.Comment.Multiline] = comment
    styles[Token.Comment.Preproc] = comment
    styles[Token.Comment.Single] = comment
    styles[Token.Comment.Special] = comment
